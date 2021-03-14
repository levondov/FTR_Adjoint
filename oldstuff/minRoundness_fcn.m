function [gg] = minRoundness(optparams)

% Compare Tom's equations to the simulation (no space charge)
%
%

% Create starting distribution
[x0,y0,xp0,yp0] = createdistribution(10000);
% x0 y0 xp0 yp0

global zv cqv sqv k_quadv k_solv k_perv
zv = []; cqv = []; sqv = []; k_quadv = []; k_solv = [];
k_perv = 0.0;

env_eqn_output_in_cart = 1;

% use the calcmoment function to calculate moments

% define variables
Q_plus_sim = 0;
Q_minus_sim = 0;
Q_x_sim = 0;
P_plus_sim = 0;
P_minus_sim = 0;
P_x_sim = 0;
E_plus_sim = 0;
E_minus_sim = 0;
E_x_sim = 0;
L_sim = 0;

% calculate all the moments for the vectors
% Q stuff
u20 = calcmoment([x0,y0],2,0);
u02 = calcmoment([x0,y0],0,2);
u11 = calcmoment([x0,y0],1,1);
Q_plus_sim = ((u20 + u02)/2.0);
Q_minus_sim = ((u20 - u02)/2.0);
Q_x_sim = (u11);
% P stuff
u11_xxp = calcmoment([x0,xp0],1,1);
u11_yyp = calcmoment([y0,yp0],1,1);
u11_yxp = calcmoment([y0,xp0],1,1);
u11_xyp = calcmoment([x0,yp0],1,1);
P_plus_sim = (u11_xxp + u11_yyp);
P_minus_sim = (u11_xxp - u11_yyp);
P_x_sim = u11_yxp + u11_xyp;
% L stuff (since we calculated what we need when doing P stuff)
L_sim = u11_xyp - u11_yxp;
% E stuff
u20 = calcmoment([xp0,yp0],2,0);
u02 = calcmoment([xp0,yp0],0,2);
u11 = calcmoment([yp0,xp0],1,1);
E_plus_sim = (u20 + u02);
E_minus_sim = (u20 - u02);
E_x_sim = 2*(u11);

y_sim = [Q_plus_sim,Q_minus_sim,Q_x_sim,P_plus_sim,P_minus_sim,P_x_sim,E_plus_sim,E_minus_sim,E_x_sim,L_sim];

% convert initial conditions from cart to larmor
%[y_lar] = cart2lar([y_sim(1,:),1.0234],-0.3073);
[y_lar] = cart2lar([y_sim,1.545],0.0);

% initial conditions for diff equations.
Q_plus = y_lar(1);
Q_minus = y_lar(2);
Q_x = y_lar(3);
P_plus = y_lar(4);
P_minus = y_lar(5);
P_x = y_lar(6);
E_plus = y_lar(7);
E_minus = y_lar(8);
E_x = y_lar(9);
L = y_lar(10);
phi_0 = y_lar(11);
%phi_0 = 1.0234;

% Initial magnet and drift positions
% [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
params = [0.50,1.0,6.4e-4,optparams(1),.054,optparams(2),optparams(3),0.054,optparams(4),optparams(5),0.054,optparams(6)];

% run integration to solve diff eqn.
h=0.001;
z_interval = [0.0,1.70]; % meters
z = z_interval(1):h:z_interval(2);
init_cond = [Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi_0]'; %Q,P,E,L
[y] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion %%
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

% transform to fixed frame
% Larmor to fixed frame transform
if env_eqn_output_in_cart
    dphi = diff(y(:,11))./diff(z');
    dphi(end+1) = dphi(end);
    y = lar2cart(y,dphi);
end


% goal is round beam.
[~,idx] = min(abs(z-params(1)));
gg = abs(y(idx,2)) + abs(y(idx,3));





end

