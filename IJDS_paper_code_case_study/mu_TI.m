function [m] = mu_TI(ws)
    % ws means wind speed
    m = 0.14*(0.75*ws+5.6)./ws;
end