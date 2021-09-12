function [z,y,motion1,motion2,params_opt] = gd_y(a,epp,emp)
%GD_F Summary of this function goes here
%   Detailed explanation goes here

% seed the random number generator and grab initial conditions
prof_offset = 0;

if nargin > 1
    [init_cond] = lse_init(epp,emp);
else
    [init_cond] = lse_init();    
end


% Initial magnet conditions
global hardedge_flag extra_params_flag extra_params_flag_firstcall
extra_params_flag = 1; % used to collect extra verbose/diagnostics data during runs.
extra_params_flag_firstcall = 1;
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    % fit parmeters for magnet length (ql) and strength (qs)
    %ql = 516.4;
    %qs = [0.001936554,0.001936481,0.001936554];
    %ql = 516.4;
    %qs = [0.001936554068875,0.001936481932222,0.001936554068875];   
    %prof_offset = 0.02157;
    ql = 100;
    qs = [.01,.01,.01];
    prof_offset = 0;    
else % for cos^2 profile
    ql = 500; %10x longer magnets
    qs = [.005,.005,.005]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end
%0.213
params_o = [0.213,0.863,-15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

params_opt = [params_o(1)*a(1),params_o(2),params_o(3)*a(2),params_o(4)*a(3),params_o(5),params_o(6)*a(4),params_o(7)*a(5),params_o(8),...
    params_o(9)*a(6),params_o(10)*a(7),params_o(11),params_o(12)*a(8),params_o(13)*a(9),params_o(14)*a(10),params_o(15)*a(11)];

params_opt'
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

if 0
    figure; 
    subplot(2,1,1); hold on;
    plot(z,EQ,'LineWidth',3); plot(z,L.^2,'LineWidth',3); plot(z,PP,'LineWidth',3);
    legend('Q_{+s} E_{+s}','L^2','P_{+s}^2');
    xlabel('Z [m]'); ylabel('a.u.');
    subplot(2,1,2); hold on;
    plot(z,motion1,'LineWidth',3); plot(z,motion2,'LineWidth',3); 
    legend('Tr[ (J_4 \Sigma)^2 ]','Det[ 2\Sigma_T ]')
    xlabel('Z [m]'); ylabel('a.u.');

    figure;
    subplot(2,1,1);
    plot(z,motion1,'LineWidth',3);
    xlabel('Z [m]'); ylabel('a.u.');
    title('0mA current: Tr[ (J_4 \Sigma)^2 ]');
    subplot(2,1,2);
    plot(z,motion2,'LineWidth',3); 
    title('0mA current: Det[ 2\Sigma_T ]')
    xlabel('Z [m]'); ylabel('a.u.');
end






