
% setup moment object
mom = MomentSolverPeriodic(10e3, 0.0e-3, [Xn_h(end,1:6)';Ep;Em;Xn_h(end,7:end)']);

% create lattice
an = ones(9,1);
mom = CreateLattice(mom, an, 1);
mom = mom.RunMoments(true);
mom.PlotBeamSize();
mom.PlotBeamEmit()
%ylim([-0.5,4.5]*1e-6);

mom = CreateLattice(mom, an, 100);
mom = mom.RunMoments(true);
mom.PlotBeamSize();
mom.PlotBeamEmit()
%ylim([-0.5,4.5]*1e-6);

%%

%X0(7) = (1/16) * ( ex2/(X0(1)+X0(2)) + ey2/(X0(1)-X0(2)) );
%X0(8) = (1/16) * ( ex2/(X0(1)+X0(2)) - ey2/(X0(1)-X0(2)) ); 




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

