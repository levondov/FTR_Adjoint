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

% setup moment object
mom = MomentSolver(10e3, 3.0e-3, X0);

% create lattice
an = ones(9,1);
mom = CreateLattice(mom, an, 1, 1);
mom.PlotBeamSize();
title('Pre optimization');
ylim([-15e-7,15e-7]);

mom = CreateLattice(mom, an_h(end,:), 1, 1);
mom.PlotBeamSize();
title('Post optimization');
ylim([-15e-7,15e-7]);


mom = CreateLattice(mom, an, 5, 1);
mom.PlotBeamSize();
title('Pre optimization');

mom = CreateLattice(mom, an_h(end,:), 5, 1);
mom.PlotBeamSize();
title('Post optimization');




%%
global extra_params_flag
extra_params_flag = 1;

an = ones(11,1); an(2)=0.0;
[z,y,m1,m2] = gd_y(an);

%%

figure; hold on;
for i = 1:10
   plot(mom.z,mom.y(:,i) - y(:,i),'.-');
end

figure; hold on;
plot(mom.y(end,:),'bo-'); plot(y(end,:),'ro');

global zv k_quadv
figure; hold on;
plot(zv,k_quadv,'.-');
plot(mom.z,mom.k,'.-');

