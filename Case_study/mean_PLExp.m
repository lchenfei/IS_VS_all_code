function [m] = mean_PLExp(ws)
    m = 2.63e-04*ws.^3-1.09e-02*ws.^2 + 1.285e-01*ws-1.32e-01;
end