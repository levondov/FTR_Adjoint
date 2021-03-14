options = optimset('Display','iter','MaxIter',250,'TolFun',1e-10);

global p1 obj1
p1 = ones(11,1); obj1 = [];

params = [0.213,0.863,15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

[O_nope,N_nope] = calcON2(z_adj,y,params,k_perv);

optparams1 = [params(1),params(3),params(4),params(6),params(7),params(9),params(10),params(12),params(13),params(14),params(15)];
optparams1 = ones(1,11);
[optparams,oo1] = fminsearch(@(params) fom_opt_fcn(params,z_adj,y_adj,y,k_perv,O_nope,N_nope),optparams1,options);

figure;
subplot(2,1,1);
plot(obj1,'linewidth',3);
subplot(2,1,2); hold on;
plot(p1'./max(max(p1)),'linewidth',2);
%%
lse_init;
% Initial magnet and drift positions
% [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
%params = [0.50,1.0,optparams(1),.01,.054,optparams(2),0.11,0.054,optparams(3),0.21,0.054, optparams(4)];
params_opt2 = [optparams(1)*params(1),params(2),optparams(2)*params(3),optparams(3)*params(4),params(5),optparams(4)*params(6),...
    optparams(5)*params(7),params(8),optparams(6)*params(9),optparams(7)*params(10),params(11),optparams(8)*params(12),...
    optparams(9)*params(13),optparams(10)*params(14),optparams(11)*params(15)];

% run integration to solve diff eqn.
rng(1234);
[init_cond] = lse_init;

global hardedge_flag k_perv

hardedge_flag = 1;
% some extra parameters we want to calculate during the simulations
global extra_params_flag extra_params_flag_firstcall
extra_params_flag = 1;
extra_params_flag_firstcall = 1;

% Initial magnet and drift positions
% [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
ql = 1;
qs = [1,1,1];
prof_offset = 0;

% run integration to solve diff eqn.
h=0.000001;
z_interval = [0.0,0.222]; % meters
z = z_interval(1):h:z_interval(2);
[y2] = ode3(@(t,Y) odefcn(t,Y,params_opt2), z_interval(1), h, z_interval(2), init_cond);
%[y1] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), -h, z_interval(1), y(end,:)');



% Constant of Motion %%
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

dphi = diff(y(:,11))./diff(z');
dphi(end+1) = dphi(end);
%y = lar2cart(y,dphi);
lseplot(z,y2,motion,'');
