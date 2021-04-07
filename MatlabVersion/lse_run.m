% Tom's envelope equations
%
% Levon D. (Oct. 2020)
%
%
%

% seed the random number generator and grab initial conditions
rng(1234);
[init_cond] = lse_init;


% Global flags for extra info
global hardedge_flag k_perv extra_params_flag extra_params_flag_firstcall
k_perv = c2perv(0); % space charge value
hardedge_flag = 1; % hard edge magent model (as opposed to cos^2 profile)
extra_params_flag = 0; % used to collect extra verbose/diagnostics data during runs.
extra_params_flag_firstcall = 0;

% Initial magnet conditions
if hardedge_flag
    % fit parmeters for magnet length (ql) and strength (qs)
    ql = 1;
    qs = [1,1,1];
    prof_offset = 0;
else % for cos^2 profile
    ql = 10; %10x longer magnets
    qs = [1,1,1]; %[0.8585,0.8363,1.1994];
    prof_offset = -0.00045;
end
params = [0.213,0.863,15e-4,... % solenoid start, length, strength
    .00425+prof_offset,.0001*ql,-18.236*qs(1),... % quad 1 start, length, strength
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),... % quad 2 start, length, strength
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),... % quad 3 start, length, strength
    45*pi/180,45*pi/180,45*pi/180]; % quad 1,2,3 rotations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% run integration to solve diff eqn.
h=0.000001; % step size
z_interval = [0.0,0.222]; % meters
z = z_interval(1):h:z_interval(2);
[y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);

% Constant of Motion %
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

dphi = diff(y(:,11))./diff(z');
dphi(end+1) = dphi(end);
%y = lar2cart(y,dphi);

lseplot(z(1:100:end),y(1:100:end,:),motion(1:100:end),'');

% reverse integration check
% [y2] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), -h, z_interval(1), y(end,:)');
% z2=z;
% lseplot(z2,flip(y2),motion,'');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% calculate adjoint related variables
dPy = zeros(3,length(z));
dEy = zeros(3,length(z));
dQy = zeros(3,length(z));
dLy = zeros(1,length(z));
FoM = zeros(1,length(z)); FoM1 = FoM; FoM2 = FoM; FoM3 = FoM; FoM4 = FoM; FoM5 = FoM;

global k_solv
k0 = 1;
komega=max(k_solv);
for i = 1:length(z) % calculate at every step, although we only need at z_end    
    % adjoint variables calculated from FoM
    dPy(:,i) = y(i,4:6);
    dEy(1,i) = k0^(-2)*(y(i,7)-0.5*komega^2*y(i,1)+k_perv)*0.5*komega^2 - 2*y(i,7)*(2*y(i,7)*y(i,1)-y(i,10)^2);
    dEy(2,i) = -k0^(2)*y(i,2);
    dEy(3,i) = -k0^(2)*y(i,3);
    dQy(1,i) = -k0^(-2)*(y(i,7)-0.5*komega^2*y(i,1)+k_perv) - 2*y(i,1)*(2*y(i,7)*y(i,1)-y(i,10)^2);
    dQy(2,i) = -k0^(-2)*y(i,8);
    dQy(3,i) = -k0^(-2)*y(i,9);
    dLy(:,i) = -2*y(i,10)*(2*y(i,7)*y(i,1)-y(i,10)^2);   
    
    % figure of merit broken into pieces for ease of reading
    FoM1(i) = sum(y(i,4:6).^2);
    FoM2(i) = (k0^2)*(y(i,2)^2 + y(i,3)^2);
    FoM3(i) = (k0^(-2))*(y(i,8)^2+y(i,9)^2);
    FoM4(i) = (k0^(-2))*(y(i,7) - 0.5*(komega^2)*y(i,1) + k_perv)^2;
    FoM5(i) = (2*y(i,7)*y(i,1) - y(i,10)^2)^2;
    FoM(i) = 0.5*(FoM1(i) + FoM2(i) + FoM3(i) + FoM4(i) + FoM5(i));
end
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

% Constant of Motion %%
L = [y_adj(:,10)];
EQ = y_adj(:,7).*y_adj(:,1) + y_adj(:,8).*y_adj(:,2) + y_adj(:,9).*y_adj(:,3);
PP = y_adj(:,4).^2 + y_adj(:,5).^2 + y_adj(:,6).^2;
motion_adj = EQ + (1/2)*L.^2 - (1/2)*PP;

lseplot(z_adj(1:100:end),y_adj(1:100:end,:),motion_adj(1:100:end),'adjoint equations');
lseplot(z_adj(1:100:end),y_ori(1:100:end,:),motion_adj(1:100:end),'original equations');

% % transform to fixed frame
% % Larmor to fixed frame transform
% if env_eqn_output_in_cart
%     dphi = diff(y(:,11))./diff(z');
%     dphi(end+1) = dphi(end);
%     y = lar2cart(y,dphi);
% end



