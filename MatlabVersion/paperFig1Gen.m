global k_perv k0 hardedge_flag
k_perv = c2perv(1.0e-3);
k0 = 10;
hardedge_flag = 1;

% seed the random number generator and grab initial conditions
Q_plus_0 = 0.5*(2.2581^2*1e-6 + 0.2258^2*1e-6);
Q_minus_0 = 0.5*(2.2581^2*1e-6 - 0.2258^2*1e-6);
Q_x_0 = 0.0;
P_plus_0 = 0;
P_minus_0 = 0;
P_x_0 = 0;
E_plus_0 = (7.0855^2*1e-6 + 0.70855^2*1e-6);
E_minus_0 = (7.0855^2*1e-6 - 0.70855^2*1e-6);
E_x_0 = 0;
L_0 = 0; 
phi_0 = 0;

init_cond = [Q_plus_0,Q_minus_0,Q_x_0,P_plus_0,P_minus_0,P_x_0,E_plus_0,E_minus_0,E_x_0,L_0,phi_0]'; %Q,P,E,L

% Initial magnet conditions
global hardedge_flag extra_params_flag extra_params_flag_firstcall
extra_params_flag = 1; % used to collect extra verbose/diagnostics data during runs.
extra_params_flag_firstcall = 1;

ql = 10;
qs = 0.1;
prof_offset = 0;    

%0.213
params_o = [0.213,0.863,-15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs,... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs,... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs,... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

params_opt = params_o;
% run integration to solve diff eqn.
h=0.00001; % step size
z_interval = [0,0.322]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode3(@(t,Y) odefcn(t,Y,params_opt), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion 1% 0.5*Tr(J_4^2 sigma^2)
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion1 = EQ + (1/2)*L.^2 - (1/2)*PP;

% Constant of Motion 2% Det(2sigma)]
motion2 = ( 2*EQ - L.^2 - PP ).^2;

%%
global zv k_quadv k_solv

figure; hold on;
plot(z,y(:,1)+y(:,2));
plot(z,y(:,1)-y(:,2));

plot(zv,k_quadv*1e-10,'k-','Linewidth',1);
plot(zv,k_solv*1e-6*0.25,'m-','Linewidth',1);

xlabel('Z position (m)'); ylabel('Moments [m^2]'); title('Beam size through FTR transformer');
legend('\langle x^2 \rangle','\langle y^2 \rangle','K_{quad}','k_{solenoid}');
grid on;