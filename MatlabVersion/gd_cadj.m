function [z,y,y_adj,params,motion] = gd_cadj(a)
%GD_CADJ Summary of this function goes here

[init_cond] = lse_init;


% Global flags for extra info
global hardedge_flag k_perv extra_params_flag extra_params_flag_firstcall

extra_params_flag = 1; % used to collect extra verbose/diagnostics data during runs.
extra_params_flag_firstcall = 1;

% Initial magnet conditions
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    %ql = 516.4;
    %qs = [0.001936554,0.001936481,0.001936554];
    ql = 100;
    qs = [.01,.01,.01];
    prof_offset = 0;
    %ql = 516.4;
    %qs = [0.001936554068875,0.001936481932222,0.001936554068875];   
    %prof_offset = 0.02157;    
else % for cos^2 profile
    ql = 500; %10x longer magnets
    qs = [.005,.005,.005]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end
%213
params_o = [0.213,0.863,-15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

params = [params_o(1)*a(1),params_o(2),params_o(3)*a(2),params_o(4)*a(3),params_o(5),params_o(6)*a(4),params_o(7)*a(5),params_o(8),...
    params_o(9)*a(6),params_o(10)*a(7),params_o(11),params_o(12)*a(8),params_o(13)*a(9),params_o(14)*a(10),params_o(15)*a(11)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% run integration to solve diff eqn.
h=0.00001; % step size
z_interval = [0.0, 0.722]; %[0.0,0.322]; % meters 0.222
z = z_interval(1):h:z_interval(2);
[y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion %
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

%y = lar2cart(y,dphi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% calculate adjoint related variables
i=length(z);
[~,~,dQy, dPy, dEy, dLy] = get_F_and_dF(y,i);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% backwards integration with adjoint variables

% run integration to solve diff eqn backwards.
idx=length(z); z_interval = [z(idx),0.0];
z_adj = z_interval(1):-h:z_interval(2);
init_cond1 = [dQy(1,idx),dQy(2,idx),dQy(3,idx),dPy(1,idx),dPy(2,idx),dPy(3,idx),dEy(1,idx),dEy(2,idx),dEy(3,idx),dLy(idx),y(idx,11)]'; %Q,P,E,L
init_cond2 = y(idx,:)';
init_cond3 = [init_cond1;init_cond2];

[y_reverse] = ode3(@(t,Y) odefcn_reverse(t,Y,params), z_interval(1), -h, z_interval(2), init_cond3 );
y_reverse1 = y_reverse(:,12:end); % original envelope backwards
y_reverse2 = y_reverse(:,1:11); % adjoint eqn.
% flip everything
y_adj = flip(y_reverse2,1);
y_ori = flip(y_reverse1,1);
z_adj = flip(z_adj);

% % transform to fixed frame
% % Larmor to fixed frame transform
% if env_eqn_output_in_cart
%     dphi = diff(y(:,11))./diff(z');
%     dphi(end+1) = dphi(end);
%     y = lar2cart(y,dphi);
% end






end

