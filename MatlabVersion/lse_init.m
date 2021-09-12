function [init_cond] = lseinit(epp,emp)

if nargin == 0
    epp = 0;
    emp = 0;
end
% Create starting distribution with 10k particles
[x0,y0,xp0,yp0] = createdistribution(10000);
% x0 y0 xp0 yp0

% extra diagnostics plotting if you want to see distribution.
if 0
figure;
subplot(2,2,1); scatter(x0,y0,'.'); ylim([-0.01,0.01]); xlim([-0.01,0.01]); title('y vs x');
subplot(2,2,2); scatter(x0,yp0,'.'); ylim([-0.01,0.01]); xlim([-0.01,0.01]); title('yp vs x');
subplot(2,2,3); scatter(y0,xp0,'.'); ylim([-0.01,0.01]); xlim([-0.01,0.01]); title('xp vs y');
subplot(2,2,4); scatter(xp0,yp0,'.'); ylim([-0.01,0.01]); xlim([-0.01,0.01]); title('yp vs xp');
end

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

% Use calculated initial conditions from Santiago insteads of calculating
% from a distribution
if 1 % using santiago values
   Q_plus_sim = 0.5*(2.2581^2*1e-6 + 0.2258^2*1e-6);
   Q_minus_sim = 0.5*(2.2581^2*1e-6 - 0.2258^2*1e-6);
   Q_x_sim = 0.0;
   P_plus_sim = 0;
   P_minus_sim = 0;
   P_x_sim = 0;
   E_plus_sim = (7.0855^2*1e-6 + 0.70855^2*1e-6);
   E_minus_sim = (7.0855^2*1e-6 - 0.70855^2*1e-6);
   %E_plus_sim = 6.0388e-05;
   %E_minus_sim = 6.9849e-05;
   E_x_sim = 0;
   L_sim = 0;
end

E_plus_sim = E_plus_sim + E_plus_sim*epp;
E_minus_sim = E_minus_sim + E_minus_sim*emp;
%E_plus_sim = E_plus_sim + E_plus_sim*-0.05;
%E_minus_sim = E_minus_sim + E_minus_sim*0.05;

y_sim = [Q_plus_sim,Q_minus_sim,Q_x_sim,P_plus_sim,P_minus_sim,P_x_sim,E_plus_sim,E_minus_sim,E_x_sim,L_sim];

% convert initial conditions from cartesian to larmor frame
[y_lar] = cart2lar([y_sim,0.0],0.0);

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
init_cond = [Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi_0]'; %Q,P,E,L


