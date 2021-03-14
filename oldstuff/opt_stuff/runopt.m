options = optimset('Display','iter','MaxIter',150,'TolFun',1e-8);
global c1
c1 = 0;
[op,oo1] = fminsearch(@(params) optgoal(params),[0.8585    0.8363    1.1994],options);

%%
rng(1234);
[init_cond] = lse_init;
ql = 10;
prof_offset = -0.00045;
params1 = [0.213,0.863,15e-4,.00425+prof_offset,.0001*ql,-18.236*op(1),0.10655+prof_offset,0.0001*ql,21.3640*op(2),0.20895+prof_offset,0.0001*ql,-18.236*op(3)];

% run integration to solve diff eqn.
h=0.00010;
z_interval = [0.0,0.212]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode4(@(t,Y) odefcn(t,Y,params1), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion %%
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

dphi = diff(y(:,11))./diff(z');
dphi(end+1) = dphi(end);
%y = lar2cart(y,dphi);
lseplot(z,y,motion);
