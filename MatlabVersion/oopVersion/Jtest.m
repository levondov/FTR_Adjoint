
verbose = false;

fixedEpm = false; % manually recalculate E each time
useJvar = true; % use the J constrain gradients
nonegative = false; % if Q+ or E+ are negative, take absolute value and continue
dontOverStepFlag = true;  % readjust gamma to make sure steps are positive

% generate some random numbers
rng('default');
rng(2);
% Initial conditions

X0 = [
   0.000046813217131
   0.000001329230425
                   0
   0.000000000000728
   0.000000000000733
                   0
   0.001276759120128
   0.000367263831509
                   0
                   0
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


%%
gammaVal = linspace(0.6,0.7,25);
gammaVal2 = 0.65;
p1 = [];
p2 = [];
for jj=1:length(gammaVal) % let it run forever and break with ctrl-c , or if conditions are met at the end of this while loop it will break automatically
    
    ii=1;
    fprintf([num2str(jj),'/',num2str(length(gammaVal))]);
    
    tmp11 = gammaVal2*df_h(:,end)';
    tmp22 = gammaVal(jj)*calcJcomponent(dj_h{end}, df_h(:,end))';
    p1(end+1) = tmp11(1);
    p2(end+1) = tmp22(1);
    
    % iterate, take a step
    tmp = Xn - tmp11 + tmp22;
    
    Xn_h(end+1,:) = tmp;
    Xval = tmp;

    mom.initialMoments = Xval';
    mom = mom.RunMoments(verbose);
    
    if useJvar
        [dj,js] = mom.CalcFoMGradientJ();
        %dj_h{end+1} = dj;
        j_h(end+1,:) = js;
    end

    % compute fom
    [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
    fprintf(['FoM: ',num2str(f_h(end)),'\n']);
    ii = ii + 1;
        
    
    % if FoM is no longer decreasing, lets recalculate the adjoint
    % equations
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    if 0
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
        if useJvar
            [dj,js] = mom.CalcFoMGradientJ();
            dj_h{end+1} = dj;
            j_h(end+1,:) = js;
        end
    end
   
end




function [dX] = calcJcomponent(djs, dw)
[~,N] = size(djs);

% calculate M matrix
M = zeros(N,N);
for i = 1:N
    for j = 1:N
        M(i,j) = dot(djs(:,i),djs(:,j));
    end
end
% inverse it
M = inv(M);

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
