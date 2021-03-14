function [params] = a2params(a)

% convert a to params

% Global flags for extra info
global hardedge_flag

% Initial magnet conditions
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    %ql = 516.4;
    %qs = [0.001936554,0.001936481,0.001936554];
    ql = 100;
    qs = [.01,.01,.01];
    prof_offset = 0;
else % for cos^2 profile
    ql = 500; %10x longer magnets
    qs = [.005,.005,.005]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end
%213
params_o = [0.213,0.863,15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

params = [params_o(1)*a(1),params_o(2),params_o(3)*a(2),params_o(4)*a(3),params_o(5),params_o(6)*a(4),params_o(7)*a(5),params_o(8),...
    params_o(9)*a(6),params_o(10)*a(7),params_o(11),params_o(12)*a(8),params_o(13)*a(9),params_o(14)*a(10),params_o(15)*a(11)];


end

