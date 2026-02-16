function [mean] = logn_mean_TI(ws,sigma)
    m = mu_TI(ws);
    mean = log(m.^2./sqrt(m.^2+sigma^2));
end