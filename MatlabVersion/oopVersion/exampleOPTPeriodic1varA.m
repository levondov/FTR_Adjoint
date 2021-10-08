
verbose = false;
plotverbose = false;
% Initial conditions
X0 = [ 
    1e-6
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0 ];

% create lattice
an = ones(5,1)';

periods = 1;

% setup moment object
mom = MomentSolverPeriodic(5e3, 0.0e-3, X0);
% create lattice
mom = CreateLattice(mom, an, periods);

% run moment + adjoint equations
mom = mom.RunMoments(verbose);
mom = mom.RunMomentsAdjoint(verbose);

% plot moments
mom.PlotBeamSize();

% Grab FoM
[f0,f0p] = mom.GetFAndDF1();

% gradient
df0 = CalcFoMGradienta(an,mom,periods);

% gradient descent parameter
gamma = (f0/sum(df0.^2));
%%
% init arrays
gamma_h = gamma;
an_h = an;
f_h = f0;
fp_h = f0p;
df_h =df0;

% adjust starting gamma
an_h(end+1,:) = an - gamma_h(end)*df0'; % iterate
mom = CreateLattice(mom, an_h(end,:), periods);
mom = mom.RunMoments(verbose);

[f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
%%
fprintf(['FoM: ',num2str(f_h(end)),'\n']);
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0;
    
    an_h(end+1,:) = an - gamma_h(end)*df0'; % iterate
    mom = CreateLattice(mom, an_h(end,:), periods);
    mom = mom.RunMoments(verbose);
    
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end
%%
if plotverbose
    figure('units','pixels','position',[100,100,1000,750]);
end
lastNpts = 100;
while 1
    ii=1;
    while f_h(end) < f_h(end-1)
        if plotverbose             
            if (length(f_h) < lastNpts)
                subplot(3,3,1:3); cla();
                plot(log10(f_h));
                for i = 4:(4+length(an)-1)
                    subplot(3,3,i); cla();
                plot(an_h(:,i-3));
                end
            else
                subplot(3,3,1:3); cla();
                plot(log10(f_h(end-100+1:end)));
                for i = 4:(4+length(an)-1)
                    subplot(3,3,i); cla();
                plot(an_h(end-100+1:end,i-3));
                end
            end
            pause(0.01);
        end
        
        fprintf(['Iterating ',num2str(ii),'\n']);        
        
        % iterate
        an_h(end+1,:) = an_h(end,:) - gamma_h(end)*df_h(:,end)';
        mom = CreateLattice(mom, an_h(end,:), periods);
        mom = mom.RunMoments(verbose);
        
        % compute fom
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        ii = ii + 1;
        
        %if (ii > 100) % if we have done more than 10 iterations without recalculating the gradient, increase gamma
        %    gamma_h(end+1) = gamma_h(end)*2.0;
        %end
    end
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    an_h(end+1,:) = an_h(end-1,:);
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    
    % calc adjoint equations
    mom = CreateLattice(mom, an_h(end,:), periods);
    mom = mom.RunMoments(verbose);
    mom = mom.RunMomentsAdjoint(verbose);
    
    % calc new gradient and gamma
    df = CalcFoMGradienta(an_h(end,:),mom, periods);
    df_h(:,end+1) = df;
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % change gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); ann = an_h(end,:);
        
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            
            an_h(end+1,:) = ann - gamma_h(end)*df_h(:,end)'; % iterate
            mom = CreateLattice(mom, an_h(end,:), periods);
            mom = mom.RunMoments(verbose);
            
            [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
            fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        end
    end
    
    %gamma = (f_h(end)/sum(df.^2))*0.01;
    % calculate learning rate
    %gamma_h(:,end+1) = gamma;
    
    if f_h(end) < 1e-20
        break;
    end
    if length(f_h) > 25000
        break;
    end
    
end





%% Helper functions

function [df] = CalcFoMGradienta(a,obj, periods)

df = zeros(length(a),1);
a_copy = a;
perturb = 0.01;

fprintf(['Calculating Gradient ']);

for i = 1:length(a)
    fprintf([num2str(i),' ']);
    
    a(i) = a_copy(i) + a_copy(i)*perturb;
    obj = CreateLattice(obj, a, periods);
    
    [O_opt,N_opt,ACT] = obj.CalcONandACTMatrix();
    for j = 1:length(obj.z)
        O_opt{j} = O_opt{j} - obj.Omat{j};
        N_opt{j} = N_opt{j} - obj.Nmat{j};
    end
    
    tmp = obj.AdjointIntegral(O_opt,N_opt);
    
    df(i) = tmp;
    
    a(i) = a_copy(i); % reset
end
fprintf('\n');

obj = CreateLattice(obj, a_copy, periods);

end

