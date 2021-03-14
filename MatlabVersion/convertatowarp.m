function [af] = convertatowarp(an)
%
%
%
%
af = zeros(15,1);

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);


% Initial magnet conditions
global hardedge_flag
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    %ql = 516.4;
    %qs = [0.001936554,0.001936481,0.001936554];
    %prof_offset = 0;
    ql = 516.4;
    qs = [0.001936554068875,0.001936481932222,0.001936554068875];   
    prof_offset = 0.02157;   
else % for cos^2 profile
    ql = 100; %10x longer magnets
    qs = [1,1,1]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end

params = [0.313,0.863,15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

af(1) = an(1)*params(1);
af(2) = params(2);
af(3) = an(2)*params(3);
af(4) = an(3)*params(4);
af(5) = params(5);
af(6) = an(4)*params(6);
af(7) = an(5)*params(7);
af(8) = params(8);
af(9) = an(6)*params(9);
af(10) = an(7)*params(10);
af(11) = params(11);
af(12) = an(8)*params(12);
af(13) = an(9)*params(13);
af(14) = an(10)*params(14);
af(15) = an(11)*params(15);



end

