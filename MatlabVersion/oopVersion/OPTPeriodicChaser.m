
flagJcorrection = false;
endOpt = false;

% generate some random numbers
rng('default');
rng(2); % use this to repeat runs, seed = 2 etc.

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

% calculate E+ and E- for a give emittance
ex2 = 6.2e-6^2; % target emittance X
ey2 = 6.2e-6^2; % target emittance Y
denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
% X0(7) = Ep;
% X0(8) = Em;
Xn = X0';

% create lattice
turns = 4;
numEvenCells = 5; % not used atm

% setup moment object
mom = MomentSolverPeriodic(10e3, 20e-3, X0); % energy, beam current, initial conditions
mom.h = 100; % integration steps in each element
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, turns, 0);

% run moment + adjoint equations
mom = mom.RunMoments();
mom = mom.RunMomentsAdjoint();

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
fbest_h = f0;
state_h = [0]; % 0 - take step, 1 - neg value gamma update, 2 - update W, 3 - update J

if 0
    % verify Jplot is correct
    bVals = linspace(-1.0,1.0,20);
    JvalsNew = [];
    for ii = 1:length(bVals)
        
        Xnew = takeStep(Xn, df0', dj0, bVals(ii));
        
        mom.initialMoments = Xnew;
        mom = mom.RunMoments(); % solve moment equations
        [dj,js] = mom.CalcFoMGradientJ();
        JvalsNew(end+1,:) = js;
        
    end
    figure; plot(bVals,JvalsNew(:,1),'.-')
end
%% adjust starting gamma
% Here we keep reepeating the same initial step and recalculating until we get a good starting gamma
% that lowers our figure of merit.

% Take a step
tmp = takeStep(Xn, df0', dj0, gamma_h(end));
Xn_h(end+1,:) = tmp;
state_h(end+1) = 0;

% Run moment equations
mom.initialMoments = Xn_h(end,:)';
mom = mom.RunMoments();

% calculate FoM
[f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
fprintf(['FoM: ',num2str(f_h(end)),'\n']);

% repeat this step until we get a FoM lower than what we started with
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0; % keep updating gamma here
    
    % Take a step
    tmp = takeStep(Xn, df0', dj0, gamma_h(end));
    Xn_h(end+1,:) = tmp;
    state_h(end+1) = 1; % update gamma and step
    
    % Run moment equations
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments();
    
    % calculate FoM
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end
fbest_h(end+1) = f_h(end);

%% Gradient descent algorithm
% Here we run our gradient descent algorithm
Xntmp = [];
while ~endOpt % let it run forever and break with ctrl-c , or if conditions are met at the end of this while loop it will break automatically
    ii=1;
    % while the FoM keeps decreasing
    while f_h(end) < f_h(end-1)
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate, take a step
        while true % while loop here is incase we get negative values
            
            % take a step
            tmp = takeStep(Xn_h(end,:), df_h(:,end)', dj_h{end}, gamma_h(end));
            
            % if the step has no negatives, break out of loop and continue
            % on
            if ( checkNegative(tmp) )
                Xn_h(end+1,:) = tmp;
                state_h(end+1) = 0;
                break;
            else
                % if step results in negative values for Q+,E+, etc, reduce
                % gamma and try again (keep looping).
                fprintf(['Neg result, attempting to reduce gamma \n']);
                gamma_h(end+1) = gamma_h(end)/2.0; % keep updating gamma here
                state_h(end+1) = 1;
                Xntmp(end+1,:) = tmp;
            end
        end
        ii = ii + 1;
        
        % correct for J deviation
        if flagJcorrection
            % iterative scheme to get J back to starting J0 value
            [tmp1, tmp2, tmp3, tmp4] = iterateJ(df_h(:,end), mom, dj_h, j_h, state_h, Xn_h);
            Xn_h = tmp1;
            dj_h = tmp2;
            j_h = tmp3;
            state_h = tmp4;
        end
        
        % Run moment equations
        mom.initialMoments = Xn_h(end,:)';
        mom = mom.RunMoments();
        
        % calculate FoM
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        
        % if we take 50+ steps without having to recompute gradients, our
        % steps are probably too small, let's start increasing gamma
        if ( ii > 50 )
            % increase gamma
            fprintf(['Too many small steps, increasing gamma... \n']);
            gamma_h(end+1) = gamma_h(end)*2.0; % keep updating gamma here
        end
        
        if (f_h(end) < fbest_h(end))
            fbest_h(end+1) = f_h(end);
        end
    end
    fprintf(['FoM increased: FoM_n: ',num2str(f_h(end)),' | FoM_n-1: ',num2str(f_h(end-1)),'\n']);
    
    % if FoM is no longer decreasing, lets recalculate the adjoint
    % equations
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    Xn_h(end+1,:) = Xn_h(end-1,:);
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    state_h(end+1) = 2;
    
    % calc adjoint equations
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments();
    mom = mom.RunMomentsAdjoint();
    
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
        jjj = 0;
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            
            % iterate, take a step
            tmp = takeStep(Xnn, df_h(:,end)', dj_h{end}, gamma_h(end));
            Xn_h(end+1,:) = tmp;
            
            % Run moment equations
            mom.initialMoments = Xn_h(end,:)';
            mom = mom.RunMoments();
            
            % calculate FoM
            [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
            fprintf(['FoM: ',num2str(f_h(end)),'\n']);
            jjj = jjj + 1;
            if (jjj == 10)
                endOpt = true;
                break;
            end
        end
        
        if (f_h(end) < fbest_h(end))
            fbest_h(end+1) = f_h(end);
        end
    end
    
    if f_h(end) < 1e-30 % break if we get this low
        endOpt = true;
    end
    if length(f_h) > 100000 % break if we have done this many iterations
        endOpt = true;
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
function [Xnh, djh, jh, sh] = iterateJ(dw, mom, djh, jh, sh, Xnh)
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
for ii = 1:5
    
    % calculate J and dJ/dX for the iteration
    mom.initialMoments = Xnh(end,:); % grab last set of moment values
    mom = mom.RunMoments(); % solve moment equations
    [dj,js] = mom.CalcFoMGradientJ(); % grab the resulting J and dJ/dW values for the given X moments
    M = calcMvector(dj, dw); % calculate M matrix
    djh{end+1} = dj; % update history
    jh(end+1,:) = js; % update history
    sh(end+1) = 3;
    
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

function [a] = calcAvector(djs, dw, b)
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
    U(i) = dot( dw, djs(:,i) );
end
end

function [M] = calcMvector(djs, dw)
[~,N] = size(djs);

% calculate M matrix
M1 = zeros(N,N);
for i = 1:N
    for j = 1:N
        M1(i,j) = dot( djs(:,i), djs(:,j) );
    end
end
% inverse it
M = inv(M1);
end

function [passFlag] = checkNegative(Xn)

passFlag = true;

if Xn(1) < 0
    passFlag = false;
end

% if (Xn(1)^2 - (Xn(2)^2 + Xn(3)^2)) < 0
%     passFlag = false;
% end

if Xn(7) < 0
    passFlag = false;
end

% if (Xn(7)^2 - (Xn(8)^2 + Xn(9)^2)) < 0
%     passFlag = false;
% end

end

