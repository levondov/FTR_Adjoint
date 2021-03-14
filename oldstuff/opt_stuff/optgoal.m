function [obj] = optgoal(inputparams)

ql = 10;
qs1=inputparams(1);
qs2=inputparams(2);
qs3=inputparams(3);
prof_offset = -0.00045;

params = [0.213,0.863,15e-4,.00425+prof_offset,.0001*ql,-18.236*qs1,0.10655+prof_offset,0.0001*ql,21.3640*qs2,0.20895+prof_offset,0.0001*ql,-18.236*qs3];

rng(1234);
[init_cond] = lse_init;

global zv cqv sqv k_quadv k_solv k_perv O N
zv = []; cqv = []; sqv = []; k_quadv = []; k_solv = [];O = {}; N = {}; k_perv = 0.0;
h=0.00010;
z_interval = [0.0,0.212]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);

obj1 = abs(y(1,1)-y(end,1)) + abs(y(1,2)-y(end,2)); + abs(y(end,3));
obj2 = 1e3*abs(y(end,10));

obj = 1e10*(obj1 + obj2);

end

