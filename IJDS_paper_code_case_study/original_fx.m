
function[p] = original_fx(X)
    wind_speed = X(:,1);
    TI = X(:,2);
    PLExp = X(:,3);
    VFlowAng = X(:,4);
    Z0 = X(:,5);
    % wind speed
    logP1 = trunc_rayleigh(wind_speed);
    % TI
    logP2 = log(lognpdf(TI,logn_mean_TI(wind_speed,0.05),logn_sd_TI(wind_speed,0.05)));
    % PLExp
    logP3 = log(normpdf(PLExp,mean_PLExp(wind_speed),sd_PLExp(wind_speed)));
    % VFlowAng
    logP4 = log(normpdf(VFlowAng,0,1) / (normcdf(10,0,1)-normcdf(-10,0,1)));
    % Z0
    logP5 = log(normpdf(Z0,0,1) / (normcdf(0.05,0,1)-normcdf(0.01,0,1)));
    p = exp(logP1+logP2+logP3+logP4+logP5);
    p(isnan(p)) = 0;
end



% density for wind speed 
function [logp] = trunc_rayleigh(x)
    Rayleigh_scale_par = sqrt(2/pi)*10;
    p = raylpdf(x,Rayleigh_scale_par)./(raylcdf(25,Rayleigh_scale_par)-raylcdf(3,Rayleigh_scale_par)).*(x>=3).*(x<=25);
    logp = log(p);
end








