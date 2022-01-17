
verbose = false;

fixedEpm = false; % manually recalculate E each time
useJvar = false; % use the J constrain gradients
nonegative = false; % if Q+ or E+ are negative, take absolute value and continue
dontOverStepFlag = false;  % readjust gamma to make sure steps are positive

% generate some random numbers
rng('default');
rng(2);
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

% create lattice
Xn = X0';

turns = 1;
numEvenCells = 5; % not used atm

% setup moment object
mom = MomentSolverPeriodic(10e3, 0.0, X0); % energy, beam current, initial conditions
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
if useJvar
    Xn_h(end+1,:) = Xn - gamma_h(end)*df0' + gamma_h(end)*calcJcomponent(dj0, df0)'; % iterate
else
    Xn_h(end+1,:) = Xn - gamma_h(end)*df0';
end
if nonegative
    Xn_h(end,1) = abs(Xn_h(end,1));
    Xn_h(end,7) = abs(Xn_h(end,7));
end
if fixedEpm
    denom1 = (Xn_h(end,1)+Xn_h(end,2)); denom2 = (Xn_h(end,1)-Xn_h(end,2));
    numer1 = ex2 + 0.5*(Xn_h(end,4)+Xn_h(end,5)).^2; numer2 = ey2 + 0.5*(Xn_h(end,4)-Xn_h(end,5)).^2;
    Ep = (numer1 / denom1) + (numer2 / denom2);
    Em = (numer1 / denom1) - (numer2 / denom2);
    Xn_h(end,7) = Ep;
    Xn_h(end,8) = Em;
end
mom.initialMoments = Xn_h(end,:)';
mom = mom.RunMoments(verbose);

% calculate FoM
[f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
fprintf(['FoM: ',num2str(f_h(end)),'\n']);

% repeat the step until we get a FoM lower than what we started with
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0; % keep updating gamma here
    
    % Take a step
    if useJvar
        Xn_h(end+1,:) = Xn - gamma_h(end)*df0' + gamma_h(end)*calcJcomponent(dj0, df0)'; % iterate
    else
        Xn_h(end+1,:) = Xn - gamma_h(end)*df0';
    end
    if nonegative
        Xn_h(end,1) = abs(Xn_h(end,1));
        Xn_h(end,7) = abs(Xn_h(end,7));
    end
    if fixedEpm
        denom1 = (Xn_h(end,1)+Xn_h(end,2)); denom2 = (Xn_h(end,1)-Xn_h(end,2));
        numer1 = ex2 + 0.5*(Xn_h(end,4)+Xn_h(end,5)).^2; numer2 = ey2 + 0.5*(Xn_h(end,4)-Xn_h(end,5)).^2;
        Ep = (numer1 / denom1) + (numer2 / denom2);
        Em = (numer1 / denom1) - (numer2 / denom2);
        Xn_h(end,7) = Ep;
        Xn_h(end,8) = Em;
    end
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments(verbose);
    
    % calculate FoM
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
end
%% Gradient descent algorithm

while 1 % let it run forever and break with ctrl-c , or if conditions are met at the end of this while loop it will break automatically
    ii=1;
    % while the FoM keeps decreasing
    while f_h(end) < f_h(end-1)    
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate, take a step
        while 1
            if useJvar
                tmp = Xn_h(end,:) - gamma_h(end)*df_h(:,end)' + gamma_h(end)*calcJcomponent(dj_h{end}, df_h(:,end))';
            else
                tmp = Xn_h(end,:) - gamma_h(end)*df_h(:,end)';
            end
            
            if dontOverStepFlag
                if ( tmp(1) < 0 || tmp(7) < 0 )
                    gamma_h(end+1) = gamma_h(end)/2.0; % gradient step too much, stepped into negative
                else
                    Xn_h(end+1,:) = tmp;
                    break;
                end
            else
                Xn_h(end+1,:) = tmp;
                break
            end
        end
        
        if nonegative
            Xn_h(end,1) = abs(Xn_h(end,1));
            Xn_h(end,7) = abs(Xn_h(end,7));
        end
        if fixedEpm
            denom1 = (Xn_h(end,1)+Xn_h(end,2)); denom2 = (Xn_h(end,1)-Xn_h(end,2));
            numer1 = ex2 + 0.5*(Xn_h(end,4)+Xn_h(end,5)).^2; numer2 = ey2 + 0.5*(Xn_h(end,4)-Xn_h(end,5)).^2;
            Ep = (numer1 / denom1) + (numer2 / denom2);
            Em = (numer1 / denom1) - (numer2 / denom2);
            Xn_h(end,7) = Ep;
            Xn_h(end,8) = Em;
        end
        mom.initialMoments = Xn_h(end,:)';
        mom = mom.RunMoments(verbose);
        if useJvar
            [dj,js] = mom.CalcFoMGradientJ();
            dj_h{end+1} = dj;
            j_h(end+1,:) = js;
        end
        
        % compute fom
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        ii = ii + 1;
        
        %if (ii > 100) % if we have done more than 10 iterations without recalculating the gradient, increase gamma
        %    gamma_h(end+1) = gamma_h(end)*2.0;
        %end
    end
    
    % if FoM is no longer decreasing, lets recalculate the adjoint
    % equations
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    Xn_h(end+1,:) = Xn_h(end-1,:);
    if nonegative
        Xn_h(end,1) = abs(Xn_h(end,1));
        Xn_h(end,7) = abs(Xn_h(end,7));
    end
    if fixedEpm
        denom1 = (Xn_h(end,1)+Xn_h(end,2)); denom2 = (Xn_h(end,1)-Xn_h(end,2));
        numer1 = ex2 + 0.5*(Xn_h(end,4)+Xn_h(end,5)).^2; numer2 = ey2 + 0.5*(Xn_h(end,4)-Xn_h(end,5)).^2;
        Ep = (numer1 / denom1) + (numer2 / denom2);
        Em = (numer1 / denom1) - (numer2 / denom2);
        Xn_h(end,7) = Ep;
        Xn_h(end,8) = Em;
    end
    f_h(end+1) = f_h(end-1);
    fp_h(end+1,:) = fp_h(end-1,:);
    
    % calc adjoint equations
    mom.initialMoments = Xn_h(end,:)';
    mom = mom.RunMoments(verbose);
    mom = mom.RunMomentsAdjoint(verbose);
    
    % calc new gradient and gamma
    df = mom.CalcFoMGradientX();
    df_h(:,end+1) = df;
    if useJvar
        [dj,js] = mom.CalcFoMGradientJ();
        dj_h{end+1} = dj;
        j_h(end+1,:) = js;
    end
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % recalculating gives no improvement, lets try adjusting gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); Xnn = Xn_h(end,:);
        
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            if useJvar
                Xn_h(end+1,:) = Xnn - gamma_h(end)*df_h(:,end)' + gamma_h(end)*calcJcomponent(dj_h{end}, df_h(:,end))'; % iterate
            else
                Xn_h(end+1,:) = Xn_h(end,:) - gamma_h(end)*df_h(:,end)';
            end
            if nonegative
                Xn_h(end,1) = abs(Xn_h(end,1));
                Xn_h(end,7) = abs(Xn_h(end,7));
            end
            if fixedEpm
                denom1 = (Xn_h(end,1)+Xn_h(end,2)); denom2 = (Xn_h(end,1)-Xn_h(end,2));
                numer1 = ex2 + 0.5*(Xn_h(end,4)+Xn_h(end,5)).^2; numer2 = ey2 + 0.5*(Xn_h(end,4)-Xn_h(end,5)).^2;
                Ep = (numer1 / denom1) + (numer2 / denom2);
                Em = (numer1 / denom1) - (numer2 / denom2);
                Xn_h(end,7) = Ep;
                Xn_h(end,8) = Em;
            end
            mom.initialMoments = Xn_h(end,:)';
            mom = mom.RunMoments(verbose);
            
            [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1(); % get FoM
            fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        end
    end
    
    %gamma = (f_h(end)/sum(df.^2))*0.01;
    % calculate learning rate
    %gamma_h(:,end+1) = gamma;
    
    if f_h(end) < 1e-30 % break if we get this low
        break;
    end
    if length(f_h) > 100000 % break if we have done this many iterations
        break;
    end    
end






function [dX] = calcJcomponent(djs, dw)
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

tmpipiece = zeros(11,1);
for i = 1:N % for each Ji
    tmpjpiece = 0;
    for j = 1:N % for each Jj
        tmpjpiece = tmpjpiece + M(i,j) * ( dot(djs(:,j),dw) );
    end
    tmpipiece = tmpipiece + djs(:,i) * tmpjpiece;
end

dX = tmpipiece;

end

