"""
The :mod:'mfis.multivariate_independent_distribution' combines a series of
indepedent continuous variable distributions to be used as an
input distribution.

@author:    D. Austin Cole <david.a.cole@nasa.gov>
            James E. Warner <james.e.warner@nasa.gov>
"""
import numpy as np
from mfis.input_distribution import InputDistribution
import scipy.stats as ss

class MVDependentDistribution(InputDistribution):
    """
    A multivariate independent distribution consisting of a series of
    independent continuous distributions. It is used to describe the
    distribution of inputs.

    Parameters
    ----------
    distributions : list
        A series of continuous distribution instances from the scipy.stats
        module. Each marginal distribution much have .rvs and .pdf functions.
    """
    def __init__(self, distributions):
        self.distributions_ = distributions


    def draw_samples(self, n_samples):
        """
        Draws and combines random input samples from the separate
        continuous distributions.

        Parameters
        ----------
        n_samples : int
            The number of samples to draw.

        Returns
        -------
        samples : array
            An n_samples by d (number of distributions) array of sample
            inputs from the Multivariate Independent distribution.
        """
        samples = np.zeros((n_samples, 5))
        print(samples)

        for i in range(n_samples):
            x_i = np.random.normal(0, 1, 3)
            x2 = np.random.normal(x_i[0], 1, 1)
            x3 = np.random.normal(x_i[0], 1, 1)
            sample_x = [x_i[0], x2[0], x3[0], x_i[1], x_i[2]]
            sample_x = np.array((sample_x))
        #print(sample_x)
            samples[i] = sample_x

        return samples

    def evaluate_pdf(self, samples):
        """
        Evaluates the probability density function of Multivariate
        Indepdendent distribution.

        Parameters
        ----------
        samples : array
            An n_samples by d array of sample inputs.

        Returns
        -------
        densities : array
            The probability densities of each sample from the Multivariate
            Independent distribution's pdf.

        """
        densities = np.ones((samples.shape[0],))

        for i in range(samples.shape[0]):
            density = ss.norm.pdf(samples[i, 0], loc = 0, scale = 1)
            density = density * ss.norm.pdf(samples[i, 1], loc = samples[i, 0], scale = 1)
            density = density * ss.norm.pdf(samples[i, 2], loc = samples[i, 0], scale = 1)
            density = density * ss.norm.pdf(samples[i, 3], loc = 0, scale = 1)
            density = density * ss.norm.pdf(samples[i, 4], loc = 0, scale = 1)
            densities[i] = density

        return densities
"""
    def ppf(self, q):
        inv_cdf = np.zeros((len(q), len(self.distributions_)))

        for i in range(len(self.distributions_)):
            inv_cdf[:, i] = self.distributions_[i].ppf(q[:,i])

        return inv_cdf 
"""
