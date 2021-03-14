function [cq,sq,Q_delta,ab4,ca_ab4,sa_ab4] = calcON_sc(z_full,Y,params,k_perv)
%CALCON Summary of this function goes here
%   Detailed explanation goes here

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);

cq = zeros(length(z_full),1);
sq = cq;
Q_delta = cq;
ab4 = cq;
ca_ab4 = cq;
sa_ab4 = cq;

% iterate through each step and calculate O(z)
for i = 1:length(z_full)
    z = z_full(i);
    
    % quads
    psi = 0;
    z_q1_start = params(4);
    z_q1_end = params(4)+params(5);
    z_q2_start = params(7);
    z_q2_end = params(7)+params(8);
    z_q3_start = params(10);
    z_q3_end = params(10)+params(11);
    
    if z >= z_q1_start && z <= z_q1_end % if inside quad 
        psi = params(9);
    end
    if z >= z_q2_start && z <= z_q2_end % if inside quad 2
        psi = params(10);
    end
    if z >= z_q3_start && z <= z_q3_end % if inside quad 3 (== quad 1)
        psi = params(11);
    end
    
    cq(i) = cos(2*Y(i,11)-2*psi);
    sq(i) = sin(2*Y(i,11)-2*psi);
    
    % space charge stuff
    Q_delta(i) = sqrt( Y(i,1)^2 - Y(i,2)^2 - Y(i,3)^2 );
    ab4(i) = 1 / Q_delta; % the 4/ab term in equation
    ca_ab4(i) = -Y(i,2) / ( (Y(i,1)+Q_delta)*Q_delta ); % 4c_alpha/ab
    sa_ab4(i) = -Y(i,3) / ( (Y(i,1)+Q_delta)*Q_delta ); % 4s_alpha/ab        
    
end

end

