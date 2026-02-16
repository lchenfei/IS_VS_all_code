#!/usr/bin/env python
# coding: utf-8

# In[3]:


"""
The :mod:'contour_adaptive_design.EntropyContourLocatorGP' conducts
adaptive design targeting a contour with the entropy
contour locator (ECL).

@author:    D. Austin Cole <david.a.cole@nasa.gov>
"""

from inspect import isfunction
import numpy as np
from pyDOE import lhs
from scipy.optimize import minimize
import scipy.stats as ss
import sys
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from tqdm import tqdm
import warnings

class ContourAdaptiveDesignGP:
    def __init__(self, gp, limit_state):
        self._initial_gp = gp
        self._alpha = gp.alpha
        self.kernel_ = gp.kernel_
        self._set_hyperparameters_from_kernel(self.kernel_)
        self.limit_state_ = limit_state
        self._K_X = self._build_covariance_matrix(gp.X_train_)
        self._Ki_X = np.linalg.inv(self._K_X)
        self.X_ = gp.X_train_
        self.y_ = np.array(gp.y_train_).reshape((-1, 1))

    def add_observation(self, x_new, y_new, update_kernel=True):
        """
        Adds new observation(s) to GP, updating relevant matrices for
        prediction.

        Parameters
        ----------
        x_new : ndarray(n,dim)
            New design location(s) (one per row).
        y_new : ndarray(n,)
            New outputs.
        update_kernel : bool, optional
            If true, the kernel is updated. The default is True.

        Returns
        -------
        None.

        """
        if update_kernel:
            self._build_new_gp(x_new, y_new)
        else:
            self._update_pred_matrices_and_data(x_new, y_new)


    def predict(self, x_pred, return_var=False):
        """
        Produces predictions from the GP's predictive equations

        Parameters
        ----------
        x_pred : ndarray(n, dim)
            A set of predictive locations, one per row.
        return_var : bool, optional
            If true, the predictive variance is also returned.
            The default is False.

        Returns
        -------
        mean: ndarray(n,)
            A set of predictive means.
        variance: ndarray(n,)
            If return_var=True, a set of predictive variances.

        """
        if x_pred.ndim == 1:
            x_pred = x_pred.reshape((1, -1))
        K_xpred_X = self.kernel_(x_pred, self.X_)
        KiY = self._Ki_X.dot(self.y_.reshape((-1, 1)))

        mean = np.dot(K_xpred_X, KiY)

        if return_var:
            var = self._calc_predictive_variance(K_xpred_X)

            return mean.flatten(), var

        return mean.flatten()


    def _build_covariance_matrix(self, X):
        K_X = self.kernel_(X)
        np.fill_diagonal(K_X, K_X.diagonal() + self._initial_gp.alpha)

        return K_X


    def _build_new_gp(self, x_new, y_new):
        new_gp = GPR(kernel=self._initial_gp.kernel, alpha=self._alpha)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            new_gp.fit(np.vstack((self.X_, x_new)), np.append(self.y_, y_new))

        self.kernel_ = new_gp.kernel_
        self._set_hyperparameters_from_kernel(new_gp.kernel_)
        self._K_X = self._build_covariance_matrix(new_gp.X_train_)
        self._Ki_X = np.linalg.inv(self._K_X)
        self.X_ = new_gp.X_train_
        self.y_ = np.array(new_gp.y_train_).reshape((-1, 1))


    def _calc_predictive_variance(self, K_xpred_X, KX_inv=None):
        if KX_inv is None:
            KX_inv = self._Ki_X
        Ki_K_xpred = KX_inv.dot(K_xpred_X.T)
        var = self._constant_value  +  self._noise_level - \
                (K_xpred_X * Ki_K_xpred.T).sum(-1)

        negative_var = (var <= 0)
        if np.sum(negative_var) > 0:
            warnings.warn(
                UserWarning("Predictive variances smaller than 0."
                            " Setting those variances to 0."))
        var[negative_var] = 0

        return var


    def _predict_variance_with_xnew(self, X_pred, x_new, X=None, KX_inv=None):
        if X is None:
            X = self.X_

        K_Xpred_Xn1 = self.kernel_(X_pred, np.vstack((X, x_new)))
        KX_inv = self._update_precision_matrix(x_new, X, KX_inv)

        K_n1i_K_Xpred = KX_inv.dot(K_Xpred_Xn1.T)
        # Diagonal of covariance matrix
        var = self._constant_value * (1 + self._alpha) + \
            self._noise_level - \
            (K_Xpred_Xn1 * K_n1i_K_Xpred.T).sum(-1)

        return var


    def _set_hyperparameters_from_kernel(self, kernel):
        kernel_params = kernel.get_params()
        self._lengthscale = \
            [value for key, value in kernel_params.items() \
             if 'length_scale' in key][0]
        self._constant_value = \
            [value for key, value in kernel_params.items() \
             if 'constant_value' in key][0]
        if any(key.endswith('noise_level') for key in kernel_params):
            self._noise_level = \
                [value for key, value in kernel_params.items() \
                 if 'noise_level' in key][0]
        else:
            self._noise_level = 0


    def _update_covariance_matrix(self, x_new, X):
        K_X_xnew = self.kernel_(X, x_new)
        k_xnew = self._build_covariance_matrix(x_new)
        np.fill_diagonal(k_xnew, k_xnew.diagonal() + self._initial_gp.alpha)

        K_Xn1 = np.append(self._K_X, K_X_xnew, axis=1)
        K_Xn1 = np.append(K_Xn1, np.hstack((K_X_xnew.T, k_xnew)), axis=0)

        return K_Xn1

    def _update_pred_matrices_and_data(self, x_new, y_new):
        self._K_X = self._update_covariance_matrix(x_new, self.X_)
        if x_new.shape[0] == 1:
            self._Ki_X = self._update_precision_matrix(x_new)
        else:
            self._Ki_X = np.linalg.inv(self._K_X)
        self.X_ = np.vstack((self.X_, x_new))
        self.y_ = np.append(self.y_, y_new)


    def _update_precision_matrix(self, x_new, X=None, K_inv=None):

        x_new = x_new.reshape((1, -1))
        if X is None:
            X = self.X_
        if K_inv is None:
            K_inv = self._Ki_X

        K_X_xnew = self.kernel_(X, x_new)
        Ki_K_xnew = K_inv.dot(K_X_xnew)

        sig2_xnew = self._constant_value* (1+self._alpha) + \
            self._noise_level - \
            np.dot(K_X_xnew.T, Ki_K_xnew)
        g_xnew = -1/sig2_xnew * Ki_K_xnew

        Ki1 = np.hstack((K_inv + g_xnew.dot(g_xnew.T)*sig2_xnew, g_xnew))
        Ki2 = np.hstack((g_xnew.T, 1/sig2_xnew))
        Kn1i = np.vstack((Ki1, Ki2))

        return Kn1i
    
    
class EntropyContourLocatorGP(ContourAdaptiveDesignGP):
    """
    Creates a GP that allows for ECL adaptive design for a targeted contour.

    Parameters
    ----------
    gp : a trained GP from the sklearn.gaussian process class
    GaussianProcessRegressor

    limit_state: the contour value to be targeted or a function
    that transforms observed responses. In the case that a function is
    supplied, the zero contour on the function's output is used.
    """
    def __init__(self, gp, limit_state):
        super().__init__(gp, limit_state)


    def calc_entropy(self, x, y_var=None, limit_state=None):
        """
        Calculates the Entropy Contour Locator values for each sample in x
        (one per row).

        Parameters
        ----------
        x : ndarray(n_samples, dim)
            An array of input samples used to evaluate ECL
        y_var : ndarray(n_samples,), optional
            An array of predictive variance estimates for the set of input
            samples *x*. This is primarily helpful for calculating lookahead
            entropy for batch selection. The default is None.
        limit_state : float or function, optional
            The contour value to be targeted or a function
            that transforms observed responses. In the case that a function is
            supplied, the zero contour on the function's output is used.
            The default is None, relying on the class attribute
            self.limit_state_.

        Returns
        -------
        entropy : ndarray(n_samples)
            An array of ECL values.

        """
        x = x.reshape(-1, self.X_.shape[1])
        if limit_state is None:
            limit_state = self.limit_state_

        if y_var is None:
            yhat, y_var = self.predict(x, return_var=True)
        else:
            yhat = self.predict(x)
        ysd = np.sqrt(y_var)

        ysd_nonzero = np.array(np.where(ysd != 0)).ravel()
        entropy = np.zeros((x.shape[0]))

        ## Only calculates ECL when predicted standard deviation is nonzero
        if ysd_nonzero.size > 0:
            if isfunction(limit_state):
                limit_state = limit_state(yhat[ysd_nonzero])
                p_less = ss.norm.cdf((limit_state)/ysd[ysd_nonzero])
            else:
                p_less = ss.norm.cdf((limit_state - yhat)/ysd[ysd_nonzero])

            p_more = 1 - p_less
            nonzero_indices, p_less_nonzero, p_more_nonzero = \
                self._find_nonzero_probabilities(p_less, p_more, ysd_nonzero)

            entropy[nonzero_indices] = \
                -p_less_nonzero*np.log(p_less_nonzero) -\
                    p_more_nonzero*np.log(p_more_nonzero)

        if x.shape[0] == 1:
            entropy = float(entropy) ## for gradient-based optimization

        return entropy


    def fit(self, n_steps, high_fidelity_model, bounds, batch_size=1,
            n_cand=None, X_cand=None, fixed_cands=False, lookahead=True,
            local_opt=True):
        """
        Performs ECL adaptive design n times, updating hyperparameters at
        each step.

        Parameters
        ----------
        n_steps : int
            Number of steps in the adaptive design.
        high_fidelity_model : class instance
            A class object that includes a predict call.
        bounds : seq
            Sequence of (min, max) pairs of the design space used as bounds
            for optimization.
        batch_size : int, optional
            Size of each batch of samples. The default is 1.
        n_cand : int, optional
            The number of candidate points to use for optimization. If None
            and X_cand is supplied, the number of rows in X_cand is used.
            Otherwise if None, 10*d is used. The default is None.
        X_cand : ndarray(n_cand, dim), optional
             A set of candidate points used for discrete optimization.
             If None, a Latin Hypercube design is used. The default is None.
        fixed_cands : bool, optional
             If True, a single candidate set is used for selecting each
             point. The default is False.
        lookahead : bool, optional
             If True, the predictive variance is updated sequentially before
             selecting the next batch point. If False, all batch samples are
             selected simultaneously. The default is True.
        local_opt : bool, optional
             If True, gradient-based optimization is performed after
             discrete optimization with X_cand. If False, only discrete
             optimization with X_cand is performed. The default is True.

        Returns
        -------
        None.

        """
        Xdim = self.X_.shape[1]
        self._check_user_inputs_for_select_samples(bounds,
                                                   n_steps*batch_size,
                                                   n_cand, X_cand,
                                                   fixed_cands)

        if X_cand is None:
            if n_cand is None:
                n_cand = 10*Xdim
            X_cand = self._generate_cand_set(n_cand, Xdim, bounds)
        else:
            n_cand = X_cand.shape[0]

        ## Performs adaptive designs steps
        for i in tqdm(range(n_steps)):
            if fixed_cands:
                if i == 0:
                    self._X_cand = X_cand
                else:
                    X_cand = self._X_cand
            else:
                if i > 0:
                    X_cand = self._generate_cand_set(n_cand, Xdim, bounds)


            ## Selects batch with discrete/continous optimization
            x_batch = self.select_samples(bounds=bounds,
                                          batch_size=batch_size,
                                          n_cand=n_cand,
                                          X_cand=X_cand,
                                          fixed_cands=fixed_cands,
                                          lookahead=lookahead,
                                          local_opt=local_opt,
                                          return_cand_preopt=False)
            y_batch = high_fidelity_model.predict(x_batch)

            self.add_observation(x_batch, y_batch, update_kernel=True)


    def opt_entropy(self, X_cand, bounds, X=None, KX_inv=None,
                    limit_state=None):
        """
        Performs gradient-based optimization ('L-BFGS-B') on ECL from a set of
        starting locations 'X_cand'. Returns the best input for a new sample.

        Parameters
        ----------
        X_cand : ndarray
            An array of starting locations for gradient-based optimization.
        bounds : seq
            Sequence of (min, max) pairs of the design space used as bounds
            for optimization.
        X : ndarray, optional
            Input locations used to calculate the covariance matrix.
            Primarily helpful for calculating look-ahead ECL.
            The default is None.
        KX_inv : TYPE, optional
            The inverse of the covariance matrix for X. Primarily helpful for
            calculating look-ahead ECL. The default is None.
        limit_state : float or function, optional
            The contour value to be targeted or a function
            that transforms observed responses. In the case that a function is
            supplied, the zero contour is used. The default is None, relying
            on the class attribute self.limit_state_.

        Returns
        -------
        x_new : ndarray(1, dim)
            Selected input location for optimal ECL.
        ECL: float
            Entropy associated with x_new

        """
        x_maybes = np.ones(X_cand.shape)
        entropy_multi = np.ones((X_cand.shape[0],))

        ## Gradient based optimization starting at each candidate point
        for i in range(X_cand.shape[0]):
            x_new = X_cand[i, :]

            res = minimize(self._neg_entropy, x_new,
                           method='L-BFGS-B', bounds=bounds,
                           args=(X, KX_inv, limit_state))
            x_maybes[i, :] = res.x
            entropy_multi[i] = res.fun

        x_new = x_maybes[np.argmin(entropy_multi), :].reshape((1, -1))

        return x_new, -np.min(entropy_multi)


    def select_samples(self, bounds, batch_size=1, n_cand=None, X_cand=None,
                       fixed_cands=False, lookahead=True, local_opt=True,
                       return_cand_preopt=False):
        """
        Selects a batch of samples using ECL adaptive design.

        Parameters
        ----------
        bounds : seq
            Sequence of (min, max) pairs of the design space used as bounds
               for optimization.
        batch_size : int, optional
            Number of samples to select. The default is 1.
        n_cand : int, optional
            The number of candidate points to use for optimization. If None
            and X_cand is supplied, the number of rows in X_cand is used.
            Otherwise if None, 10*d is used. The default is None.
        X_cand : ndarray(n_cand, dim), optional
            A set of candidate points used for discrete optimization. If None,
            a Latin Hypercube design is used. The default is None.
        fixed_cands : bool, optional
            If True, a single candidate set is used for selecting each
            point in the batch. The default is False.
        lookahead : bool, optional
            If True, the predictive variance is updated sequentially before
            selecting the next batch point. If False, all batch samples are
            selected simultaneously. The default is True.
        local_opt : bool, optional
            If True, gradient-based optimization is performed after discrete
            optimization with X_cand. If False, only discrete
            optimization with X_cand is performed. The default is True.
        return_cand_preopt : bool, optional
            If True, the selected candidate point in each discrete
            optimization is returned in addition to the final optimized
            samples. This is only pertinent when 'local_opt=True'.
            The default is False.

        Returns
        -------
        opt_samples: ndarray(batch_size, dim)
            Batch of samples selected using ECL optimization.
            If local_opt=True, the samples are based on a two-stage
            optimization: discrete followed by gradient-based.
            If local_opt=False, samples are based on discrete optimization.
        cand_samples: ndarray(batch_size, dim)
            If return_cand_preopt=True and local_opt=True, a batch of
            candidate points selected during the discrete optimization stage.

        """

        X_cand, n_cand = \
            self._check_user_inputs_for_select_samples(bounds, batch_size,
                                                       n_cand, X_cand,
                                                       fixed_cands)

        if not local_opt and return_cand_preopt:
            warnings.warn('Local optimization not performed. '
                          'Both sets of samples are equivalent.')

        if lookahead and batch_size > 1:
            samples = \
                self._select_batch_with_lookahead(bounds=bounds,
                                                  n_cand=n_cand,
                                                  X_cand=X_cand,
                                                  batch_size=batch_size,
                                                  fixed_cands=fixed_cands,
                                                  local_opt=local_opt,
                                                  return_cand_preopt=\
                                                      return_cand_preopt)
        else:
            samples = \
                self._select_batch_without_lookahead(bounds=bounds,
                                                     n_cand=n_cand,
                                                     X_cand=X_cand,
                                                     batch_size=batch_size,
                                                     fixed_cands=fixed_cands,
                                                     local_opt=local_opt,
                                                     return_cand_preopt=\
                                                         return_cand_preopt)

        return samples


    def _check_user_inputs_for_select_samples(self, bounds, n_select, n_cand,
                                              X_cand, fixed_cands):
        Xdim = self.X_.shape[1]
        if X_cand is not None:
            if X_cand.shape[1] != Xdim:
                print('The number of columns in X_cand and X do not match.')
                sys.exit(0)

        if len(bounds) != Xdim:
            print('The length of bounds is less than the number of columns '
                  'in X.')
            sys.exit(0)

        if X_cand is None:
            if n_cand is None:
                n_cand = 10*Xdim
            X_cand = self._generate_cand_set(n_cand, Xdim, bounds)

        else:
            n_cand = X_cand.shape[0]

        if fixed_cands and n_cand < n_select:
            print('The number of candidate points is less than the '
                  'batch size.')
            sys.exit(0)

        return X_cand, n_cand



    def _find_nonzero_probabilities(self, p_less, p_more, ysd_nonzero):
        calc_indices = np.logical_and(p_less > 1e-300, p_more > 1e-300)
        calc_indices = np.where((p_less > 1e-300) & (p_more > 1e-300))

        nonzero_indices = np.array(ysd_nonzero).ravel()[calc_indices]
        p_less_sub = p_less[calc_indices]
        p_more_sub = p_more[calc_indices]

        return nonzero_indices, p_less_sub, p_more_sub


    def _find_opt_bounds(self, xstar, X, bounds):

        bound_box = np.zeros((X.shape[1], 2))
        ## Determines bounding box for local optimization based on 'buffer'
        ## closest points in X
        buffer = 2
        for i in range(X.shape[1]):
            sortedXi = X[X[:, i].argsort(), i]
            ii = np.asarray(np.where(sortedXi == xstar[0, i]))
            reached_max = ii + buffer - (len(X) - 1)
            if reached_max > 0:
                bound_box[i, 1] = bounds[i][1]
            else:
                bound_box[i, 1] = sortedXi[ii + buffer]

            reached_max = ii - buffer
            if reached_max < 0:
                bound_box[i, 0] = bounds[i][0]
            else:
                bound_box[i, 0] = sortedXi[ii - buffer]

        return tuple(map(tuple, bound_box))


    def _generate_cand_set(self, n, dim, bounds):
        X_cand = lhs(dim, n)
        for k in range(dim):
            X_cand[:, k] = (bounds[k][1] - bounds[k][0])*X_cand[:, k] + \
                bounds[k][0]

        return X_cand


    def _neg_entropy(self, x_pred, X=None, KX_inv=None, limit_state=None):
        if X is None:
            X = self.X_
        x_pred = x_pred.reshape((-1, X.shape[1]))

        if KX_inv is None:
            KX_inv = self._Ki_X

        if X is None and KX_inv is None:
            K_xpred_X = self.kernel_(x_pred, X)
            y_pred_var = self._calc_predictive_variance(K_xpred_X, KX_inv)
        else:
            y_pred_var = None

        entropy = self.calc_entropy(x_pred, y_pred_var, limit_state)

        return -entropy

    def _remove_near_replicate(self, x_opt, x_cand, X_batch):
        dists = np.linalg.norm(X_batch-x_opt, axis=1)
        if np.min(dists) < 1e-4:
            return x_cand
        return x_opt

    def _select_batch_without_lookahead(self, bounds, n_cand, X_cand,
                                        batch_size, fixed_cands=False,
                                        local_opt=True,
                                        return_cand_preopt=False):
        xdim = self.X_.shape[1]
        if fixed_cands or batch_size == 1:
            X_batch, cand_index = \
                self._select_n_samples(batch_size, bounds, X_cand,
                                       local_opt=local_opt)
            X_batch_cand = X_cand[cand_index, :].reshape(batch_size, xdim)
            self._X_cand = np.delete(X_cand, cand_index, axis=0)

        else:
            X_batch = np.zeros((batch_size, xdim))
            X_batch_cand = np.zeros((batch_size, xdim))

            ## Select remaining points in batch
            for i in range(batch_size):
                if i > 0:
                    X_cand = self._generate_cand_set(n_cand, xdim, bounds)

                x_new, cand_index = \
                    self._select_n_samples(n=1, bounds=bounds, X_cand=X_cand,
                                           local_opt=local_opt)
                X_batch_cand[i, :] = X_cand[cand_index, :]

                if local_opt:
                    X_batch[i, :] = \
                        self._remove_near_replicate(x_new, X_batch_cand[i, :],
                                                    X_batch)
                else:
                    X_batch[i, :] = X_batch_cand[i, :]


        if return_cand_preopt:
            return X_batch, X_batch_cand
        else:
            return X_batch



    def _select_batch_with_lookahead(self, bounds, n_cand, X_cand,
                                     batch_size=5, fixed_cands=False,
                                     local_opt=True, return_cand_preopt=False):

        X = np.copy(self.X_)
        X_batch = np.zeros((batch_size, X.shape[1]))
        X_batch_cand = np.zeros((batch_size, X.shape[1]))
        KX_inv = self._Ki_X

        ## Select points in batch
        for i in range(batch_size):
            if not fixed_cands and i > 0:
                X_cand = self._generate_cand_set(n_cand, X.shape[1], bounds)

            x_new, cand_index = \
                self._select_n_samples(n=1, bounds=bounds, X_cand=X_cand,
                                       X=X, KX_inv=KX_inv, local_opt=local_opt)
            X_batch_cand[i, :] = X_cand[cand_index, :]

            if local_opt:
                X_batch[i, :] = \
                    self._remove_near_replicate(x_new, X_batch_cand[i, :],
                                                X_batch)
            else:
                X_batch[i, :] = X_batch_cand[i, :]

            if fixed_cands:
                X_cand = np.delete(X_cand, cand_index, axis=0)

            if i != batch_size-1:
                x_new_opt = (X_batch[i-1, :]).reshape((1, -1))
                KX_inv = self._update_precision_matrix(x_new_opt, X, KX_inv)
                X = np.vstack((X, x_new_opt))

        if fixed_cands:
            self._X_cand = X_cand

        if return_cand_preopt:
            return X_batch, X_batch_cand

        return X_batch


    def _select_candidate_samples(self, n, X_cand, X=None, KX_inv=None,
                                  limit_state=None):
        """
        Selects the n best samples through discrete optimization of ECL
        with a single candidate set

        Parameters
        ----------
        n : int
            The number of samples to select.
        X_cand : ndarray(n_cand, dim)
            A candidate set of input locations from which to select new inputs.
        X : ndarray, optional
            Input locations used to calculate the covariance matrix.
            Primarily helpful for calculating look-ahead ECL.
            The default is None.
        KX_inv : TYPE, optional
            The inverse of the covariance matrix for X. Primarily helpful for
            calculating look-ahead ECL. The default is None.
        limit_state : float or function, optional
            The contour value to be targeted or a function
            that transforms observed responses. In the case that a function is
            supplied, the zero contour is used. The default is None, relying
            on the class attribute self.limit_state_.

        Returns
        -------
        x_new : ndarray(n, dim)
            Set of input samples selected from candidate points.

        """
        ## Calculate y_var for lookahead ECL
        if len(X_cand) > n:
            if KX_inv is not None:
                K_xpred_X = self.kernel_(X_cand, X)
                y_var = self._calc_predictive_variance(K_xpred_X, KX_inv)
            else:
                y_var = None

            entropy_step_i = self.calc_entropy(X_cand, y_var, limit_state)

            if np.all(entropy_step_i == 0):
                warnings.warn("All candidate entropy values are zero.")
                ind = np.arange(len(entropy_step_i))
                np.random.shuffle(ind)
                selected_ind = ind[0:n]
            else:
                selected_ind = (-entropy_step_i).argsort()[:n]
        else:
            selected_ind = np.arange(n)

        return selected_ind

    def _select_opt_samples(self, n, X_start, bounds, X=None, KX_inv=None):
        x_new_opt = np.empty(X_start.shape)
        ecl_values = np.empty((X_start.shape[0],))
        for i in range(n):
            x_start_i = X_start[i, :].reshape((1, -1))
            x_new_opt[i, :], ecl_values[i] = \
                self.opt_entropy(x_start_i, bounds, X=X, KX_inv=KX_inv)

        selected_ind = (-ecl_values).argsort()[:n]
        x_new = x_new_opt[selected_ind, :].reshape(n, X_start.shape[1])

        return x_new

    def _select_n_samples(self, n, bounds, X_cand, X=None, KX_inv=None,
                          local_opt=True):

        """
        Selects n new input samples by first selecting the n best candidate
        points and then performing gradient-based optimization beginning
        at each selected candidate.

        Parameters
        ----------
        n : int
            The number of samples to select.
        X_cand : ndarray(n_cand, dim)
            A candidate set of input locations from which to select new inputs.
        X : ndarray, optional
            Input locations used to calculate the covariance matrix.
            Primarily helpful for calculating look-ahead ECL.
            The default is None.
        bounds : seq
            Sequence of (min, max) pairs of the design space used as bounds
            for optimization.
        KX_inv : TYPE, optional
            The inverse of the covariance matrix for X. Primarily helpful for
            calculating look-ahead ECL. The default is None.
        limit_state : float or function, optional
            The contour value to be targeted or a function
            that transforms observed responses. In the case that a function is
            supplied, the zero contour is used. The default is None, relying
            on the class attribute self.limit_state_.
        return_cand_preopt : bool, optional
            If true, the candidate points from the discrete optimization steps
            are returned. The default is False.

        Returns
        -------
        xnew: ndarray(n, dim)
            Selected input samples.

        xnew_cand: ndarray(n, dim)
            Selected candidate points from discrete optimization.

        """

        selected_cand_ind = self._select_candidate_samples(n, X_cand, X,
                                                           KX_inv)
        x_new_cand = X_cand[selected_cand_ind, :].reshape(n, self.X_.shape[1])

        if local_opt:

            x_new_opt = self._select_opt_samples(n, x_new_cand, bounds, X=X,
                                                 KX_inv=KX_inv)
            if n > 1:
                for i in range(1, n):

                    x_new_opt[i, :] = self._remove_near_replicate(
                        x_new_opt[i, :], x_new_cand[i, :],
                        (x_new_opt[0:i, :]).reshape(i, self.X_.shape[1]))

            return x_new_opt, selected_cand_ind

        return x_new_cand, selected_cand_ind


    
#############
# Example
###########
class BraninHooModel:
    def __init__(self):
        pass

    def predict(self, X):
        X0 = np.copy(X)
        
        a = 1
        # b = 5.1/(4*np.pi**2)
        b = 5/(4*np.pi**2)
        c = 5/np.pi
        r = 6
        s = 10
        t = 1/(8*np.pi)
        
        # y = a*(X0[:,1]-b*X0[:,0]**2+c*X0[:,0]-r)**2 + \
        #     s*(1-t)*np.sin(X0[:,0]) + s
        y = a*(X0[:,1]-b*X0[:,0]**2+c*X0[:,0]-r)**2 + \
            s*(1-t)*np.cos(X0[:,0]) + s
        
        return y

class HartmannModel:
    def __init__(self):
        self.C_ =  [1, 1.2, 3, 3.2]
 
        self.a_ = np.array([[10, .05, 3, 17],[3, 10, 3.5, 8],
                            [17, 17, 1.7, .05], [3.5, .1, 10, 10], 
                            [7, 8, 17, .1], [8, 14, 8, 14]])
        
        self.p_ = np.array([[.1312, .2329, .2348, .4047], 
                            [.1696, .4135, .1451, .8828],
                            [.5569, .8307, .3522, .8732], 
                            [.0124, .3736, .2883, .5743],
                            [.8283, .1004, .3047, .1091], 
                            [.5886, .9991, .6650, .0381]])


    def predict(self, X):
        y = np.zeros((X.shape[0],))
        for i in range(4):
            expon = np.zeros((X.shape[0],))
            for j in range(6):
                expon += self.a_[j,i]*(X[:,j] - self.p_[j,i])**2
            y += self.C_[i]*np.exp(-expon)

        return y
    

class IshigamiModel:
    def __init__(self, a_coef=5, b_coef=0.1):
        self.a_coef = a_coef
        self.b_coef = b_coef


    def predict(self, X):
        outputs = np.sin(X[:, 0]) + self.a_coef*np.sin(X[:, 1])**2 + \
            self.b_coef*X[:, 2]**4 * np.sin(X[:, 0])

        return outputs

    
class Numerical_Ex1:
    def __init__(self):
        pass

    def predict(self, X):
        X0 = np.copy(X)
        mu = 40*(1-np.exp(-0.2*np.sqrt((X0[:,0]**2+X0[:,1]**2)/2))) \
            + 20*(1-np.exp(-0.2*np.abs(X[:,0]))) + 5*(1-np.exp(-0.2*np.sqrt((X0[:,1]**2+X0[:,2]**2+X0[:,3]**2)/3))) \
            - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,2])) \
            - np.exp(np.cos(2*np.pi*X0[:,1]*X0[:,2])) \
            - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,1]))
        y = mu + np.random.normal(0,1,X0.shape[0])
        #(Pt,"0.01" = 18.99,"0.05" = 15.02, "0.1" = 12.88)
        
        return y

class Numerical_Ex2:
    def __init__(self):
        pass

    def predict(self, X):
        X0 = np.copy(X)
        mu = 40*(1-np.exp(-0.2*np.sqrt((X0[:,0]**2+X0[:,1]**2)/2))) \
        + 20*(1-np.exp(-0.2*np.abs(X0[:,0])))+ 5*(1-np.exp(-0.2*np.sqrt((X0[:,1]**2+X0[:,2]**2)/2))) \
        - 0.1*np.exp(-0.2 * np.sqrt((X0[:,3]**2+X0[:,4]**2)/2)) \
        - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,2])) \
        - np.exp(np.cos(2*np.pi*X0[:,1]*X0[:,2])) \
        - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,1]))
        y = mu + np.random.normal(0,1,X0.shape[0])
        #(Pt,"0.01" = 18.74,"0.05" = 14.79, "0.1" = 12.63)
        
        return y
    
class Numerical_Ex3:
    def __init__(self):
        pass

    def predict(self, X):
        X0 = np.copy(X)
        mu = 40*(1-np.exp(-0.2*np.sqrt((X0[:,0]**2+X0[:,1]**2)/2))) \
        + 20*(1-np.exp(-0.2*np.abs(X0[:,0])))+ 5*(1-np.exp(-0.2*np.sqrt((X0[:,1]**2+X0[:,2]**2)/2))) \
        - 0.1*np.exp(-0.2 * np.sqrt((X0[:,3]**2+X0[:,4]**2)/2)) - 0.1*np.exp(-0.2 * np.sqrt((X0[:,5]**2+X0[:,6]**2)/2)) \
        - 0.1*np.exp(-0.2 * np.sqrt((X0[:,7]**2+X0[:,8]**2+X0[:,9]**2)/3)) \
        - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,2])) \
        - np.exp(np.cos(2*np.pi*X0[:,1]*X0[:,2])) \
        - np.exp(np.cos(2*np.pi*X0[:,0]*X0[:,1]))
        y = mu + np.random.normal(0,1,X0.shape[0])
        #(Pt,"0.01" = 18.74,"0.05" = 14.79, "0.1" = 12.63)
        
        return y

