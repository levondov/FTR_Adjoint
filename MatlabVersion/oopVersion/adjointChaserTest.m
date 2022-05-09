% Test all the functions of adjoint chaser

% Initial conditions
% 0mA
X0 = [
   0.000407212723699
  -0.000183484865165
                   0
   0.000000005174747
   0.000000006065963
                   0
   0.118445742459199
   0.053370142482708
                   0
                   0
                   0] * 1e-3;  
               
Y0 = X0;               

%% Create moment objects for each X & Y

momX = MomentSolverPeriodicChaser(10e3, 0, X0); momX.h = 100; % default is h=10000 which is too much for now
momY = MomentSolverPeriodicChaser(10e3, 0, Y0); momY.h = 100;

%% Create our lattice of 16 periods 

% here we have 32 quads each with a adjustable strength dB, thus our a vector is length 32
period = 16;
qlength = 0.02;
qstart = [0.0, 0.03, 0.07];
db = [0.5, -0.5, 0.5];
qend = [qstart(1)+qlength/2.0, qstart(2)+qlength, qstart(3)+qlength/2.0];
qrot = [0.0, 0.0, 0.0, 0.0];
zstart = 0.0;
zend = 0.08;

% setup lattice
momX = momX.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);
momY = momY.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

%% Run Moment equations

% run moment + adjoint equations
momX = momX.RunMoments();
momX = momX.RunMomentsAdjoint();
momY = momY.RunMoments();
momY = momY.RunMomentsAdjoint();

%% Gradients of W
% we want to calculate dW/dX and dW/da both for X and Y

% dW/dX and dW/dY
dWdX = momX.GetWGradientX(); % returns dW/dX for Q,P,E,L,phi
dWdY = momY.GetWGradientX(); % note function name has an "X", but this is for our Y variable

%dW(X)/da and dW(Y)/da
dWdaX = momX.GetWGradientA(); % heavy operation, takes ~9 seconds each
dWdaY = momY.GetWGradientA();

%% Run Moment equations for FE

% run moments + adjoint equations FE
momX = momX.RunMoments();
momX = momX.RunMomentsAdjointFE();
momY = momY.RunMoments();
momY = momY.RunMomentsAdjointFE();

%% Gradients of FE
% we want to calculate dFE/dX and dFE/da both for X and Y

% dW/dX and dW/dY
dFEdX = momX.GetFEGradientX(); % returns dFE/dX for Q,P,E,L,phi
dFEdY = momY.GetFEGradientX(); % note function name has an "X", but this is for our Y variable

%dW(X)/da and dW(Y)/da
dFEdaX = momX.GetFEGradientA(); % heavy operation, takes ~9 seconds each
dFEdaY = momY.GetFEGradientA();



