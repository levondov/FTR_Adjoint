options = optimset('Display','iter','MaxIter',50);

global p1 obj1
p1 = zeros(11,1); obj1 = [];
params = [0.213,0.863,15e-4,.00425+prof_offset,.0001*ql,-18.236*qs(1),0.10655+prof_offset,0.0001*ql,...
    21.3640*qs(2),0.20895+prof_offset,0.0001*ql,-18.236*qs(3), 45*pi/180, 45*pi/180, 45*pi/180];

optparams1 = [params(1),params(3),params(4),params(6),params(7),params(9),params(10),params(12),params(13),params(14),params(15)];
[optparams,oo1] = fminsearch(@sc_opt_fcn,optparams1,options);

figure; 
subplot(2,1,1);
plot(obj1,'linewidth',3);
subplot(2,1,2); hold on;
plot(p1'./max(max(p1)),'linewidth',2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization , only run once
%
% initial seed and distributions
rng(1234);
[init_cond] = lse_init;

% setup flags for the solver
global hardedge_flag k_perv

% some extra parameters we want to calculate during the simulations
global extra_params_flag extra_params_flag_firstcall
extra_params_flag = 1;
extra_params_flag_firstcall = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% run integration to solve diff eqn.
h=0.00001;
z_interval = [0.0,0.212]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode4(@(t,Y) odefcn(t,Y,optparams), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion %%
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

dphi = diff(y(:,11))./diff(z');
dphi(end+1) = dphi(end);
%y = lar2cart(y,dphi);
lseplot(z,y,motion,'rotated');
