%rre = (rand(1,11)*2-1)*0.25;
%%
%
%
%
% Gradient descent with adjoint stuff
global k_perv k0 hardedge_flag
k_perv = c2perv(1e-3);
k0 = 10;
hardedge_flag = 1;
flag_noangle = 0;

%an = [1.40, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

an = ones(1,11);
an(1) = 1.2;

% compute adjoint equations
[z,y,y_adj,params,motion] = gd_cadj(an);
[O_nope,N_nope] = calcON2(z,y,params,k_perv);

% compute initial conditions for FoM and dFoM
[f0,f0p] = gd_F(an); f00 = f0;
df0 = gd_dF(an,z,y,y_adj,k_perv,O_nope,N_nope);
ss = 1;
gamma = (f0/sum(df0.^2));

%%
gamma_h = gamma;
an_h = an; 
f_h = f0; 
fp_h = f0p;
df_h =df0; 
ss_h = ss;

% adjust starting gamma
an_h(end+1,:) = an - ss*gamma_h(end)*df0'; % iterate
if flag_noangle
    an_h(end,9:11) = 1; %% ANGLE
end
[f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:)); % get FoM
fprintf(['FoM 0: ',num2str(f_h(end-1)),'\n']);
fprintf(['FoM 1: ',num2str(f_h(end)),'\n']);

%%

while 1
    if f_h(end) < (f_h(end-1)/2.0) % case 1       
        % do nothing
        fprintf(['F1 < F0/2 , no gamma change','\n']);
    else % case 2, f1>f0/2
        fprintf(['F1 > F0/2 , yes gamma change','\n']);
        loop_flag = 1;
        sm = 1; fm = f_h(end); fm0 = f_h(end-1); fmp = fp_h(end,:); anm = an_h(end-1,:);
        while loop_flag            
            % compute new step size
            sm = (0.5*(sm^2)*fm0)/(fm - fm0*(1-sm));

            % compute new FoM with stepsize
            an_new = an_h(end-1,:) - sm*gamma_h(end)*df_h(:,end)'; % update
            if flag_noangle
                an_new(1,9:11) = 1; %% ANGLE
            end        
            [fm_new,fp_new] = gd_F(an_new); % get FoM 
            fprintf(['Trying s_m+1: ',num2str(sm),' | FoM_m+1: ',num2str(fm_new),'\n']);

            % go through cases
            if fm0 > fm_new
                if fm_new > fm
                   % pick am and fm
                   an_h(end+1,:) = anm;
                   f_h(end+1) = fm; fp_h(end+1,:) = fmp;
                   fprintf(['Accepted m | FoM_m: ',num2str(fm),'\n']);
                else
                    % pick am_new and fm_new
                    an_h(end+1,:) = an_new;
                    f_h(end+1) = fm_new; fp_h(end+1,:) = fp_new;
                    fprintf(['Accepted m+1 | FoM_m+1: ',num2str(fm_new),'\n']);
                end
                loop_flag = 0; % stop looping.
            else
                % increment step size and try again
                fm = fm_new;
                fmp = fp_new;
                anm = an_new;
            end        
        end
        ss_h(end+1) = sm; % take updated step size        
    end
    
    fprintf(['\n']);
    fprintf(['Recalculating adjoint eqn ... \n']);
     % recompute gradient
    [z,y,y_adj,params] = gd_cadj(an_h(end,:));
    [O_nope,N_nope] = calcON2(z,y,params,k_perv);
    df = gd_dF(an_h(end,:),z,y,y_adj,k_perv,O_nope,N_nope);
    df_h(:,end+1) = df;

    % iterate
    gamma_h(end+1) = f_h(end)/sum(df.^2);
    an_h(end+1,:) = an_h(end,:) - ss*gamma_h(end)*df';
    if flag_noangle
        an_h(end,9:11) = 1; %% ANGLE
    end

    % calculate new FoM
    [f_h(end+1),fp_h(end+1,:)] = gd_F(an_h(end,:)); % get FoM 
    fprintf(['FoM 0: ',num2str(f_h(end-1)),'\n']);
    fprintf(['FoM 1: ',num2str(f_h(end)),'\n']);
    
end

%%






