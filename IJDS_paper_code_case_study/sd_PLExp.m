function [sd] = sd_PLExp(ws)
    sd = 7.767e-05*ws.^3 - 3.43e-03*ws.^2 + 3.4e-02*ws + 1.3e-01;
end