
flagJcorrection = true;

% generate some random numbers
rng('default');
rng(2); % use this to repeat runs, seed = 2 etc.

pamount = 10.0;
% You can adjust the number of points for the grid scan so it takes less
% time.
NNp = linspace(-0.5,0.5,40); % perturb Q+
NNm = linspace(-0.5,0.5,40); % perturb Q-
FoMresults = zeros(length(NNp),length(NNm));

        
% create lattice
turns = 36; % number of periods/turns

% solution we are using
X0 = [
   0.007363052819400
  -0.003144679253067
                   0
   0.000000302378449
   0.000000350230487
                   0
   0.554383041347526
   0.239127669593029
                   0
                   0
                   0]' * 1e-4;       
               
% setup moment object
mom = MomentSolverPeriodic(10e3, 20e-3, X0); % energy, beam current, initial conditions
mom.h = 100; % integration steps in each element
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, turns, 0);

% grid scan
for i = 1:length(NNp)
    for j = 1:length(NNm)
        % update initial conditions
        Xn = X0;
        Xn(1) = X0(1)*(1+NNp(i));
        Xn(2) = X0(2)*(1+NNm(j));
        
        % run moment + adjoint equations           
        mom.initialMoments = Xn;
        mom = mom.RunMoments();
                
        % Grab FoM
        [f0,f0p] = mom.GetFAndDF1();
        FoMresults(i,j) = f0;
        
        fprintf([num2str(i),'/',num2str(length(NNp)),' | ',num2str(j),'/',num2str(length(NNm)),'\n']);

    end
end

%% Plot the grid

figure;
pcolor(NNp,NNm,log10(FoMresults));
xlabel('Q+ perturb'); ylabel('Q- perturb');
colorbar();
colormap(jet);
shading interp;
title('FoM vs perturbations 20mA');
savefig('test');

%% Plot moments for check

% create lattice
turns = 128;
numEvenCells = 5; % not used atm
pqp = -0.27;
pqm = -0.25;

X0 = [
   0.007363052819400*(1+pqp)
  -0.003144679253067*(1+pqm)
                   0
   0.000000302378449
   0.000000350230487
                   0
   0.554383041347526
   0.239127669593029
                   0
                   0
                   0]' * 1e-4;                
% setup moment object
mom = MomentSolverPeriodic(10e3, 20e-3, X0); % energy, beam current, initial conditions
mom.h = 100; % integration steps in each element
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, turns, 0);
mom = mom.RunMoments();
mom.PlotBeamSize();

