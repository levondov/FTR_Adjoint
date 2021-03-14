function [k_perv] = c2perv(current)
% convert from beam current to beam perveance value

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);
v = beta*c;

k_perv = (1/(4*pi))*(c*377*current) / (m*v^3*gamma^3/e);

end

%%% 0.021289 * current