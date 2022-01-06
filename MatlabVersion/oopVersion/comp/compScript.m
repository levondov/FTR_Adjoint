% data file info
% run1 - single period warp solution
% run2 - 20 periods

%% Load WARP data
fname = 'data/run2.mat';
dt = load(fname);
dt_warp = zeros(length(dt.data{1}),10);
dt_warp_rms = zeros(length(dt.data{1}),3);
for ii = 1:10
   dt_warp(:,ii) = dt.data{3+ii}'; 
end
for ii = 1:3
   dt_warp_rms(:,ii) = dt.data{ii}';     
end
% get rid of trailing zeros
if ( find(dt_warp(:,1) == 0,1) )
    idx = find(dt_warp(:,1) == 0,1) - 1;
    tmp = dt_warp(1:idx,:);
    tmp2 = dt_warp_rms(1:idx,:);
    dt_warp = tmp;
    dt_warp_rms = tmp2;
end

%% python moment data

fname = 'data/run1py.mat';
dt = load(fname);
dt_py = dt.data{2}';
dt_py_rms = dt.data{1}';


%% matlab moment data
% Initial conditions
X0 = dt_warp(1,:); X0(end+1) = 0.0;
% setup moment object
mom = MomentSolverPeriodic(10e3, 0.0e-3, X0);
mom.h = 0.001;

% create lattice
periods = 20;
mom = CreateLattice(mom, an, periods);
mom = mom.RunMoments(true);
dt_matlab = zeros(length(mom.z),10);
dt_matlab_rms = zeros(length(mom.z),3);
for ii = 1:10
   dt_matlab(:,ii) = mom.y(:,ii); 
end
dt_matlab_rms(:,1) = mom.z;
dt_matlab_rms(:,2) = mom.y(:,1) + mom.y(:,2);
dt_matlab_rms(:,3) = mom.y(:,1) - mom.y(:,2);

%% Compare
plotEveryX = 1;

figure('units','pixels','position',[200,300,1600,400]);

subplot(1,4,1); hold on;
plot(dt_py_rms(:,1),dt_py(:,1),'k-');
plot(dt_py_rms(:,1),dt_py(:,2),'r-');
plot(dt_py_rms(:,1),dt_py(:,3),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,1),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,2),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,3),'b*--');
legend('Q+ matlab','Q- matlab','Qx matlab','Q+ warp','Q- warp','Qx warp','NumColumns',2,'Location','southoutside');
title('Q');


subplot(1,4,2); hold on;
plot(dt_py_rms(:,1),dt_py(:,4),'k-');
plot(dt_py_rms(:,1),dt_py(:,5),'r-');
plot(dt_py_rms(:,1),dt_py(:,6),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,4),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,5),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,6),'b*--');
legend('P+ matlab','P- matlab','Px matlab','P+ warp','P- warp','Px warp','NumColumns',2,'Location','southoutside');
title('P');


subplot(1,4,3); hold on;
plot(dt_py_rms(:,1),dt_py(:,7),'k-');
plot(dt_py_rms(:,1),dt_py(:,8),'r-');
plot(dt_py_rms(:,1),dt_py(:,9),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,7),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,8),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,9),'b*--');
legend('E+ matlab','E- matlab','Ex matlab','E+ warp','E- warp','Ex warp','NumColumns',2,'Location','southoutside');
title('E');


subplot(1,4,4); hold on;
plot(dt_matlab_rms(:,1),dt_matlab(:,10),'k-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,10),'k*--');
legend('L matlab','L warp','NumColumns',1,'Location','southoutside');
title('L');

stophere();
%% Compare
plotEveryX = 3;

figure('units','pixels','position',[200,300,1600,400]);

subplot(1,4,1); hold on;
plot(dt_matlab_rms(:,1),dt_matlab(:,1),'k-');
plot(dt_matlab_rms(:,1),dt_matlab(:,2),'r-');
plot(dt_matlab_rms(:,1),dt_matlab(:,3),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,1),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,2),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,3),'b*--');
legend('Q+ matlab','Q- matlab','Qx matlab','Q+ warp','Q- warp','Qx warp','NumColumns',2,'Location','southoutside');
title('Q');


subplot(1,4,2); hold on;
plot(dt_matlab_rms(:,1),dt_matlab(:,4),'k-');
plot(dt_matlab_rms(:,1),dt_matlab(:,5),'r-');
plot(dt_matlab_rms(:,1),dt_matlab(:,6),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,4),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,5),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,6),'b*--');
legend('P+ matlab','P- matlab','Px matlab','P+ warp','P- warp','Px warp','NumColumns',2,'Location','southoutside');
title('P');


subplot(1,4,3); hold on;
plot(dt_matlab_rms(:,1),dt_matlab(:,7),'k-');
plot(dt_matlab_rms(:,1),dt_matlab(:,8),'r-');
plot(dt_matlab_rms(:,1),dt_matlab(:,9),'b-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,7),'k*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,8),'r*--');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,9),'b*--');
legend('E+ matlab','E- matlab','Ex matlab','E+ warp','E- warp','Ex warp','NumColumns',2,'Location','southoutside');
title('E');


subplot(1,4,4); hold on;
plot(dt_matlab_rms(:,1),dt_matlab(:,10),'k-');
plot(dt_warp_rms(1:plotEveryX:end,1),dt_warp(1:plotEveryX:end,10),'k*--');
legend('L matlab','L warp','NumColumns',1,'Location','southoutside');
title('L');


%% Plot errors

figure; hold on;
linestyle = '-';
for ii = 1:10
    if ii >= 7
        linestyle = '--';
    end
   plot(dt_matlab_rms(:,1),abs(dt_matlab(:,ii)-dt_warp(:,ii))./dt_matlab(1,ii), 'linewidth',2,'linestyle',linestyle);
   %plot(dt_matlab_rms(:,1),abs(dt_matlab(:,ii)-dt_warp(:,ii)), 'linewidth',2,'linestyle',linestyle); 
end
xlim([0,dt_matlab_rms(end-10,1)]);
legend('Q+','Q-','Qx','P+','P-','Px','E+','E-','Ex','L');

%%

fname = 'data/run4_5sc.mat';
% load warp stuff
dt = load(fname);
aa = 0.25*dt.data{14}.^2;
bb = 0.25*dt.data{16}.^2;
zz = 0:0.0001:1.52;

dt = load(fname);
dt_warp = zeros(length(dt.data{1}),10);
dt_warp_rms = zeros(length(dt.data{1}),3);
for ii = 1:10
   dt_warp(:,ii) = dt.data{3+ii}'; 
end
for ii = 1:3
   dt_warp_rms(:,ii) = dt.data{ii}';     
end

% matlab moment data
% Initial conditions
X0 = dt_warp(1,:); X0(end+1) = 0.0;
% setup moment object
Energy = 10e3;
gamma = 1+((Energy)/(510998.9461));
mom = MomentSolverPeriodic(Energy, 5.0e-3 / gamma^2, X0);
mom.h = 1000;

% create lattice
periods = 60;
mom = CreateLattice(mom, an, periods);
mom = mom.RunMoments(true);
aam = mom.y(:,1) + mom.y(:,2);
bbm = mom.y(:,1) - mom.y(:,2);

evx = 20;
figure; hold on;
plot(mom.z,aam,'k');
plot(zz(1:evx:end),aa(1:evx:end),'k*');

plot(mom.z,bbm,'r');
plot(zz(1:evx:end),bb(1:evx:end),'r*');







