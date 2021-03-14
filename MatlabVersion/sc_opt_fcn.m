function [int_value] = sc_opt_fcn(params)
%FOM_OPT_FCN Summary of this function goes here
%   Detailed explanation goes here
global p1 obj1
p1(:,end+1) = params';

% initial seed and distributions
rng(1234);
[init_cond] = lse_init;

% setup flags for the solver
global hardedge_flag k_perv
global y_nosc
k_perv = c2perv(1e-3); % 1 mA
hardedge_flag = 1;
% some extra parameters we want to calculate during the simulations
global extra_params_flag extra_params_flag_firstcall
extra_params_flag = 0;
extra_params_flag_firstcall = 1;

ql = 1;
qs = [1,1,1];
prof_offset = 0;

params_o = [0.213,0.863,15e-4,.00425+prof_offset,.0001*ql,-18.236*qs(1),0.10655+prof_offset,0.0001*ql,...
    21.3640*qs(2),0.20895+prof_offset,0.0001*ql,-18.236*qs(3),45*pi/180,45*pi/180,45*pi/180];

params_new = [params_o(1)*params(1),0.863,params_o(3)*params(2),params_o(4)*params(3),.0001,params_o(6)*params(4),params_o(7)*params(5),0.0001,...
    params_o(9)*params(6),params_o(10)*params(7),0.0001,params_o(12)*params(8),params_o(13)*params(9),params_o(14)*params(10),params_o(15)*params(11)];


% run integration to solve diff eqn.
h=0.000001;
z_interval = [0.0,0.222]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode4(@(t,Y) odefcn(t,Y,params_new), z_interval(1), h, z_interval(2), init_cond);

y=y(end,:);
goal = sum(abs(y(end,:)-y_nosc(end,:)));
int_value = 1e10*goal;


obj1(end+1) = int_value;
end

