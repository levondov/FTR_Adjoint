

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

% setup moment object
mom = MomentSolver(5e3, 3.0e-3, X0);


% create lattice
an = ones(9,1);
mom = CreateLattice(mom, an);

% run moment solution
mom = mom.RunMoments(true);
mom = mom.RunMomentsAdjoint(true);

mom.PlotBeamSize();







%% Helper functions

% define a lattice
function obj = CreateLattice(obj, a)

    db = [-.18236*a(2), .213640*a(4), -.18236*a(6)];
    qstart = [0.00425*a(1), 0.10655*a(3), 0.20895*a(5)];
    qend = [qstart(1)+0.01, qstart(2)+0.01, qstart(3)+0.01];
    qrot = [45.0*pi/180*a(7), 45.0*pi/180*a(8), 45.0*pi/180*a(9)];
    
    obj = obj.CreateLatticeProfile(db,qstart,qend,qrot, 0, 0.322, 1, true);
    
end

% gradient of F with respect to parameters
function [df, obj] = CalcAdjointIntegral(a,obj)

df = zeros(length(a),1);
a_copy = a;

fprintf(['Calculating Gradient ']);

for i = 1:length(a)
    fprintf([num2str(i),' ']);
    
    a(i) = a_copy(i) + a_copy(i)*perturb;
    obj = CreateLattice(obj, a);

    [O_opt,N_opt] = obj.CalcONMatrix();
    for j = 1:length(z)
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

    obj = CreateLattice(obj, a_copy);

end

