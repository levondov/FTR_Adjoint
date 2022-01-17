
% generate some random numbers
rng('default');
rng(2); % use this to repeat runs, seed = 2 etc.

% Initial conditions
X0 = [
    rand % Q+
    rand % Q-
    0    % Qx
    rand % P+
    rand % P-
    0    % Px
    0    % E+
    0    % E-
    0    % Ex
    0    % L
    0 ] * 1e-3; % phi (leftover from Flat-Round stuff, represents rotation angle in larmor frame)

% calculate E+ and E- for a give emittance
ex2 = 7.6e-6^2; % target emittance X
ey2 = 7.6e-6^2; % target emittance Y
denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
X0(7) = Ep;
X0(8) = Em;
Xn = X0';

% create lattice
turns = 1;
numEvenCells = 5; % not used atm

% setup moment object
mom = MomentSolverPeriodic(10e3, 5.0e-3, X0); % energy, beam current, initial conditions
mom.h = 1000; % integration steps in each element
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, turns, 0);

% run moment + adjoint equations
mom = mom.RunMoments(verbose);
mom = mom.RunMomentsAdjoint(verbose);

% plot moments
mom.PlotBeamSize();

% Grab FoM
[f0,f0p] = mom.GetFAndDF1();

% gradient
df0 = mom.CalcFoMGradientX();
[dj0,j0] = mom.CalcFoMGradientJ();

% gradient descent parameter
gamma = (f0/sum(df0.^2));

% init arrays
gamma_h = gamma;
Xn_h = Xn;
f_h = f0;
fp_h = f0p;
df_h = df0;
dj_h = {dj0};
j_h = j0;

%% adjust starting gamma
% Here we keep reepeating the same initial step and recalculating until we get a good starting gamma
% that lowers our figure of merit.

% Take a step
tmp = takeStep(Xn, df0', dj0, gamma_h(end));
Xn_h(end+1,:) = tmp;

% Run moment equations
mom.initialMoments = Xn_h(end,:)';
mom = mom.RunMoments(verbose);

% calculate FoM
[f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
fprintf(['FoM: ',num2str(f_h(end)),'\n']);

% repeat this step until we get a FoM lower than what we started with
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0; % keep updating gamma here
    
    % Take a step
    tmp = takeStep(Xn, df0', dj0, gamma_h(end));
    Xn_h(end+1,:) = tmp;
    
    % Run moment equations
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments(verbose);
    
    % calculate FoM
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end

%% Gradient descent algorithm
% Here we run our gradient descent algorithm

while 1 % let it run forever and break with ctrl-c , or if conditions are met at the end of this while loop it will break automatically
    ii=1;
    % while the FoM keeps decreasing
    while f_h(end) < f_h(end-1)
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate, take a step
        tmp = takeStep(Xn_h(end,:), df_h(:,end)', dj_h{end}, gamma_h(end));
        Xn_h(end+1,:) = tmp;
        
        % iterative scheme to get J back to starting J0 value
        [tmp1, tmp2, tmp3] = iterateJ(df_h(:,end), mom, dj_h, j_h, Xn_h);
        Xn_h = tmp1;
        dj_h = tmp2;
        j_h = tmp3;
        
        % Run moment equations
        mom.initialMoments = Xn_h(end,:)';
        mom = mom.RunMoments(verbose);
        
        % calculate FoM
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
    end
    
    % if FoM is no longer decreasing, lets recalculate the adjoint
    % equations
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    Xn_h(end+1,:) = Xn_h(end-1,:);
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    
    % calc adjoint equations
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments(verbose);
    mom = mom.RunMomentsAdjoint(verbose);
    
    % calc new gradient and gamma
    df = mom.CalcFoMGradientX();
    df_h(:,end+1) = df;
    
    % calc new J and gradient J
    [dj,js] = mom.CalcFoMGradientJ();
    dj_h{end+1} = dj;
    j_h(end+1,:) = js;
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % recalculating gives no improvement, lets try adjusting gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); Xnn = Xn_h(end,:);
        
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            
            % iterate, take a step
            tmp = takeStep(Xnn, df_h(:,end)', dj_h{end}, gamma_h(end));
            Xn_h(end+1,:) = tmp;
            
            % Run moment equations
            mom.initialMoments = Xn_h(end,:)';
            mom = mom.RunMoments(verbose);
            
            % calculate FoM
            [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
            fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        end
    end
    
    if f_h(end) < 1e-30 % break if we get this low
        break;
    end
    if length(f_h) > 100000 % break if we have done this many iterations
        break;
    end
end


%
% Take a gradient step
function [X] = takeStep(X0, dw, dj, b)

% grab "a" piece
a = calcAvector(dj, dw, b);

%sum J component
Jpart = 0;
for i = 1:length(a)
    Jpart = Jpart + a(i)*dj(:,i);
end

% take step
X = X0 - b*dw - Jpart';
end

%
% Iterative scheme to fix J
function [Xnh, djh, jh] = iterateJ(dw, mom, djh, jh, Xnh)
% 
% dw - gradient of W (dW/dX)
% mom - moment solver object
% djh - history list of dJ/dX 
% jh - history list of J
% Xnh - history list of X values 

% grab starting j0, the ideal j0 that we want to have.
j0 = jh(1,:);
mu = 1;
[~,N] = size(j0);

% do 10 iterations
for ii = 1:10
    
    % calculate J and dJ/dX for the iteration
    mom.initialMoments = Xnh(end,:); % grab last set of moment values
    mom = mom.RunMoments(); % solve moment equations
    [dj,js] = mom.CalcFoMGradientJ(); % grab the resulting J and dJ/dW values for the given X moments
    M = calcMvector(dj, dw); % calculate M matrix
    djh{end+1} = dj; % update history
    jh(end+1,:) = js; % update history
    
    % calculate sums for the iteration scheme
    % Xn = Xn-1 + \sum a_i,n-1 * dJ_idX (see Tom's equations)
    fullVal = 0;
    for i = 1:N
        aVal = 0;
        for j = 1:N
            aVal = aVal + M(i,j)*( j0(j)-js(j) );
        end
        fullVal = fullVal + aVal*dj(:,i);
    end
    
    % iterate
    Xnh(end+1,:) = Xnh(end,:) + mu*fullVal';
end

end

function [U] = calcAvector(djs, dw, b)
[~,N] = size(djs);
a = zeros(N,1);

U = calcUvector(djs, dw);
M = calcMvector(djs, dw);

% calculate a
for i = 1:N
    sumComp = 0;
    for j = 1:N
        sumComp = sumComp + M(i,j)*U(j);
    end
    a(i) = -b * sumComp;
end

end

function [U] = calcUvector(djs, dw)
[~,N] = size(djs);

U = zeros(N,1);

for i = 1:N
    U(i) = dot( dw,djs(:,i) );
end
end

function [M] = calcMvector(djs, dw)
[~,N] = size(djs);

% calculate M matrix
M1 = zeros(N,N);
for i = 1:N
    for j = 1:N
        M1(i,j) = dot( djs(:,i),djs(:,j) );
    end
end
% inverse it
M = inv(M1);
end

