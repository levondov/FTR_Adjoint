function [int_value] = fom_opt_fcn(params,z_adj,y_adj,y,k_perv,O_nope,N_nope)
%FOM_OPT_FCN Summary of this function goes here
%   Detailed explanation goes here
global p1 obj1
p1(:,end+1) = [params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9),params(10),params(11)]';

ql = 1;
qs = [1,1,1];
prof_offset = 0;

params_o = [0.213,0.863,15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

params_opt = [params_o(1)*params(1),0.863,params_o(3)*params(2),params_o(4)*params(3),.0001,params_o(6)*params(4),params_o(7)*params(5),0.0001,...
    params_o(9)*params(6),params_o(10)*params(7),0.0001,params_o(12)*params(8),params_o(13)*params(9),params_o(14)*params(10),params_o(15)*params(11)];

[O_opt,N_opt] = calcON2(z_adj,y,params_opt,k_perv);
for j = 1:length(z_adj)
    O_opt{j} = O_opt{j} - O_nope{j};
    N_opt{j} = N_opt{j} - N_nope{j};
end
[int_value] = find_int(z_adj,y_adj,y,O_opt,N_opt);

int_value = 1/(1e10*abs(int_value));


obj1(end+1) = int_value;
end

