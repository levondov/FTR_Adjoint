%
%
%
%
% Gradient descent with adjoint stuff
global k_perv
k_perv = c2perv(0);

an = ones(1,11);
an = ones(1,11) + rre;
%an = an*rand(11)*0.1;
% compute adjoint equations
[z,y,y_adj,params] = gd_cadj(an);
[O_nope,N_nope] = calcON2(z,y,params,k_perv);
%%
% compute initial conditions for FoM and dFoM
f0 = gd_F(an); f00 = f0;
df0 = gd_dF(an,z,y,y_adj,k_perv,O_nope,N_nope);
%gamma = (f0/sum(df0.^2))*0.1;
gamma = sqrt(sum(df0.^2))/sum(df0.^2)*0.03;


%%
gamma_h = gamma; 
an_h = an; 
f_h = f0; 
df_h =df0; 

an_h(end+1,:) = an_h(end,:) - gamma*df0'; % iterate
f_h(end+1) = gd_F(an_h(end,:)); % get FoM

while 1
    ii=1;
    while f_h(end) <= f_h(end-1)
        % grab gradient
        df = gd_dF(an_h(end,:),z,y,y_adj,k_perv,O_nope,N_nope);
        df_h(:,end+1) = df;
        
        num = ((an_h(end,:)-an_h(end-1,:))*(df_h(:,end)-df_h(:,end-1)))*0.1;
        gamma = abs(num)/sum((df_h(:,end)-df_h(:,end-1)).^2);
        gamma_h(end+1) = gamma;
        
        % iterate
        fprintf(['Iterating ',num2str(ii),'\n']);
        an_h(end+1,:) = an_h(end,:) - gamma*df';
        
        % compute fom
        f_h(end+1) = gd_F(an_h(end,:));
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        ii = ii + 1;
    end

    % recompute adjoint equation stuff for new direction
    fprintf(['Changing directions, recomputing adjoint equations \n']);
    
    an_h(end+1,:) = an_h(end-1,:);        
    f_h(end+1) = gd_F(an_h(end,:));
    
    [z,y,y_adj,params] = gd_cadj(an_h(end,:));
    [O_nope,N_nope] = calcON2(z,y,params,k_perv);
    
end


%%

an_sol = an_h(7,:);

% Initial magnet and drift positions
% [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
%params = [0.50,1.0,optparams(1),.01,.054,optparams(2),0.11,0.054,optparams(3),0.21,0.054, optparams(4)];
params_opt2 = [an_sol(1)*params(1),params(2),an_sol(2)*params(3),an_sol(3)*params(4),params(5),an_sol(4)*params(6),...
    an_sol(5)*params(7),params(8),an_sol(6)*params(9),an_sol(7)*params(10),params(11),an_sol(8)*params(12),...
    an_sol(9)*params(13),an_sol(10)*params(14),an_sol(11)*params(15)];

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