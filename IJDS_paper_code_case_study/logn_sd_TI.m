function [sd] = logn_sd_TI(ws,sigma)
    m = mu_TI(ws);
    sd = sqrt(log(1+sigma^2./m.^2));
end