%rre = (rand(1,11)*2-1)*0.25;
%%
%
%
%

% Gradient descent with adjoint stuff
global k_perv k0 hardedge_flag e1 e2
k_perv = c2perv(0.0);
k0 = 10;
e1 = 0.0;
e2 = 10.0;
hardedge_flag = 1;
flag_noangle = 0;
flag_nosoldmove = 0;

%an = [1.40, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
%1.020022195126994

an = [
   1.157924056225655
  -0.977797295545406
   0.999016701746482
   1.014754055333640
   1.021900004464416
   1.042514483002313
   1.010851419901695
   1.044696927788185
   1.030980974725414
   1.031091113146560
   1.020180786900699]';    

an = [1.15, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

if flag_noangle
    an(9:11) = 1.0;
end
if flag_nosoldmove
    an(1) = 1.0;
end

%an = ones(1,11); an(1) = 1.10;
%an([2,4,6,8]) = 0.01;

% compute adjoint equations
[z,y,y_adj,params,motion] = gd_cadj(an);
[O_nope,N_nope] = calcON2(z,y,params,k_perv);

% compute initial conditions for FoM and dFoM
[f0,f0p] = gd_F(an); f00 = f0;
df0 = gd_dF(an,z,y,y_adj,k_perv,O_nope,N_nope);
gamma = (f0/sum(df0.^2));


%%
gamma_h = gamma;
an_h = an; 
f_h = f0; 
fp_h = f0p;
df_h =df0; 

% adjust starting gamma
an_h(end+1,:) = an - an.^2 * gamma_h(end) * df0'; % iterate
if flag_noangle
    an_h(end,9:11) = 1; %% ANGLE
end
if flag_nosoldmove
    an_h(end,1) = an(1);
end
[f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:)); % get FoM
fprintf(['FoM: ',num2str(f_h(end)),'\n']);
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0;
    an_h(end+1,:) = an - gamma_h(end)*df0'; % iterate
    if flag_noangle
        an_h(end,9:11) = 1; %% ANGLE
    end
    if flag_nosoldmove
        an_h(end,1) = an(1);
    end
    [f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:)); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end
%%
while 1
    ii=1;
    while f_h(end) < f_h(end-1)        
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate
        %an_h(end+1,:) = an_h(end,:) - gamma_h(end)*df_h(:,end)';
        an_h(end+1,:) = an_h(end,:) - gamma_h(end)*df_h(:,end)';
        if flag_noangle
            an_h(end,9:11) = 1; %% ANGLE
        end
        if flag_nosoldmove
            an_h(end,1) = an(1);
        end
        
        % compute fom
        [f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:));
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        ii = ii + 1;
        
        if (ii > 20) % if we have done more than 10 iterations without recalculating the gradient, increase gamma
            gamma_h(end+1) = gamma_h(end)*2.0;
        end
    end

    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    an_h(end+1,:) = an_h(end-1,:);  
    if flag_noangle
        an_h(end,9:11) = 1; %% ANGLE
    end
    if flag_nosoldmove
        an_h(end,1) = an(1);
    end
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    
    % calc adjoint equations
    [z,y,y_adj,params] = gd_cadj(an_h(end,:));
    [O_nope,N_nope] = calcON2(z,y,params,k_perv);
    
    % calc new gradient and gamma
    df = gd_dF(an_h(end,:),z,y,y_adj,k_perv,O_nope,N_nope);
    df_h(:,end+1) = df;
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % change gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); ann = an_h(end,:);
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            an_h(end+1,:) = ann - gamma_h(end)*df_h(:,end)'; % iterate
            if flag_noangle
                an_h(end,9:11) = 1; %% ANGLE
            end
            if flag_nosoldmove
                an_h(end,1) = an(1);
            end            
            [f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:)); % get FoM
            fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        end
    end
    
    %gamma = (f_h(end)/sum(df.^2))*0.01;
    % calculate learning rate
    %gamma_h(:,end+1) = gamma;
    
    if f_h(end) < 1e-15
       break; 
    end
    if length(f_h) > 5000
        break;
    end
    
end