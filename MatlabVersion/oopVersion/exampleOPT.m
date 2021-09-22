
verbose = true;
% Initial conditions
X0 = [ 0.5*(2.2581^2*1e-6 + 0.2258^2*1e-6)
   0.5*(2.2581^2*1e-6 - 0.2258^2*1e-6)
   0.0
   0.0
   0.0
   0.0
   (7.0855^2*1e-6 + 0.70855^2*1e-6)
   (7.0855^2*1e-6 - 0.70855^2*1e-6)
   0.0
   0.0
   0.0 ];

% create lattice
an = ones(9,1)';

% setup moment object
mom = MomentSolver(5e3, 0.0e-3, X0);

% create lattice
mom = CreateLattice(mom, an);

% run moment + adjoint equations
mom = mom.RunMoments(verbose);
mom = mom.RunMomentsAdjoint(verbose);

mom.PlotBeamSize();
pause(0.01);

% Grab FoM
[f0,f0p] = mom.GetFAndDF1();

% gradient
df0 = CalcFoMGradient(an,mom);

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
mom = CreateLattice(mom, an_h(end,:));
mom = mom.RunMoments(verbose);

[f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM

fprintf(['FoM: ',num2str(f_h(end)),'\n']);
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0;
    
    an_h(end+1,:) = an - gamma_h(end)*df0'; % iterate
    mom = CreateLattice(mom, an_h(end,:));
    mom = mom.RunMoments(verbose);
    
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end

while 1
    ii=1;
    while f_h(end) < f_h(end-1)
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate
        an_h(end+1,:) = an_h(end,:) - gamma_h(end)*df_h(:,end)';
        mom = CreateLattice(mom, an_h(end,:));
        mom = mom.RunMoments(verbose);
        
        % compute fom
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
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
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    
    % calc adjoint equations
    mom = CreateLattice(mom, an_h(end,:));
    mom = mom.RunMoments(verbose);
    mom = mom.RunMomentsAdjoint(verbose);
    
    % calc new gradient and gamma
    df = CalcFoMGradient(an_h(end,:),mom);
    df_h(:,end+1) = df;
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % change gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); ann = an_h(end,:);
        
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            
            an_h(end+1,:) = ann - gamma_h(end)*df_h(:,end)'; % iterate
            mom = CreateLattice(mom, an_h(end,:));
            mom = mom.RunMoments(verbose);
            
            [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
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





%% Helper functions

% define a lattice
function obj = CreateLattice(obj, a)

    db = [-.18236*a(2), .213640*a(4), -.18236*a(6)];
    qstart = [0.00425*a(1), 0.10655*a(3), 0.20895*a(5)];
    qend = [qstart(1)+0.01, qstart(2)+0.01, qstart(3)+0.01];
    qrot = [45.0*pi/180*a(7), 45.0*pi/180*a(8), 45.0*pi/180*a(9)];
    
    obj = obj.CreateLatticeProfile(db,qstart,qend,qrot, 0, 0.322, 1, false);
    
end

% gradient of F with respect to parameters
function [df] = CalcFoMGradient(a,obj)

df = zeros(length(a),1);
a_copy = a;
perturb = 0.01;

fprintf(['Calculating Gradient ']);

for i = 1:length(a)
    fprintf([num2str(i),' ']);
    
    a(i) = a_copy(i) + a_copy(i)*perturb;
    obj = CreateLattice(obj, a);

    [O_opt,N_opt] = obj.CalcONMatrix();
    for j = 1:length(obj.z)
        O_opt{j} = O_opt{j} - obj.Omat{j};
        N_opt{j} = N_opt{j} - obj.Nmat{j};
    end
    
    tmp = obj.AdjointIntegral(O_opt,N_opt);

    %f5_tmp = ( y(ii,7) + 0.5*(komega^2)*y(ii,1) - komega*y(ii,10) );
    %f4_tmp = ( y(ii,7) - 0.5*(komega^2)*y(ii,1) + k_perv );
    %tmp5 = k0^(-2) * f5_tmp * e2 * ( 0.5*y(ii,1) - 0.5*y(ii,10)*(1/komega) );
    %tmp4 = k0^(-2) * f4_tmp * e1 * ( -0.5*y(ii,1) );
    df(i) = tmp;% + tmp5 + tmp4;% / ((k0^2) * (y(1,1)^2));
    %df(i) = tmp / (a_copy(i)*perturb);
    
    a(i) = a_copy(i); % reset
end
    fprintf('\n');

    obj = CreateLattice(obj, a_copy);

end

