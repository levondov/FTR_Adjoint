%%

% Initial conditions
X0 = [
   0.001252188797105
  -0.000682443064303
                   0
   0.000000000264310
   0.000000000088867
                   0
   0.131234365905479
   0.071522747222643
                   0
                   0
                   0 ] * 1e-3;

% calculate E+ and E-
ex2 = 7.6e-6^2;
ey2 = 7.6e-6^2;
denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
Ep = (numer1 / denom1) + (numer2 / denom2);
Em = (numer1 / denom1) - (numer2 / denom2);
X0(7) = Ep;
X0(8) = Em;

% create lattice
Xn = X0';

periods = 1;

% setup moment object
mom = MomentSolverPeriodic(10e3, 5e-3, X0);
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, periods);

%% vary step size and see results

% change step size
mom.h = 0.0001;

% run moment
mom = mom.RunMoments(verbose);

% look at moments and see how much they change for different step size
mom.y % contains all 11 moments as a function of z through the lattice