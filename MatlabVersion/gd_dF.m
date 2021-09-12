function [df] = gd_dF(a,z,y,y_adj,k_perv,O_nope,N_nope)
%GD_DF Summary of this function goes here
%   Detailed explanation goes here

perturb = 0.01;
% Initial magnet conditions
global hardedge_flag k0
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    % fit parmeters for magnet length (ql) and strength (qs)
    %ql = 516.4;
    %qs = [0.001936554,0.001936481,0.001936554];
    ql = 100;
    qs = [.01,.01,.01];
    prof_offset = 0;
    %ql = 516.4;
    %qs = [0.001936554068875,0.001936481932222,0.001936554068875];   
    %prof_offset = 0.02157;    
else % for cos^2 profile
    ql = 500; %10x longer magnets
    qs = [.005,.005,.005]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end

params_o = [0.213,0.863,-15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

df = zeros(length(a),1);
a_copy = a;
global k0 k_solvv norm_flag e1 e2 k_perv
global k_solvv
komega = k_solvv;
fprintf(['Calculating Gradient ']);
ii=length(z);
for i = 1:length(a)
    fprintf([num2str(i),' ']);
    a(i) = a_copy(i) + a_copy(i)*perturb;
    params_opt = [params_o(1)*a(1),params_o(2),params_o(3)*a(2),params_o(4)*a(3),params_o(5),params_o(6)*a(4),params_o(7)*a(5),params_o(8),...
        params_o(9)*a(6),params_o(10)*a(7),params_o(11),params_o(12)*a(8),params_o(13)*a(9),params_o(14)*a(10),params_o(15)*a(11)];    
    
    [O_opt,N_opt] = calcON2(z,y,params_opt,k_perv);
    for j = 1:length(z)
        O_opt{j} = O_opt{j} - O_nope{j};
        N_opt{j} = N_opt{j} - N_nope{j};
    end
    
    [tmp] = find_int(z,y_adj,y,O_opt,N_opt);

    f5_tmp = ( y(ii,7) + 0.5*(komega^2)*y(ii,1) - komega*y(ii,10) );
    f4_tmp = ( y(ii,7) - 0.5*(komega^2)*y(ii,1) + k_perv );
    tmp5 = k0^(-2) * f5_tmp * e2 * ( 0.5*y(ii,1) - 0.5*y(ii,10)*(1/komega) );
    tmp4 = k0^(-2) * f4_tmp * e1 * ( -0.5*y(ii,1) );
    df(i) = tmp + tmp5 + tmp4;% / ((k0^2) * (y(1,1)^2));
    %df(i) = tmp / (a_copy(i)*perturb);
    
    a(i) = a_copy(i); % reset
end
fprintf(['\n']);



end

