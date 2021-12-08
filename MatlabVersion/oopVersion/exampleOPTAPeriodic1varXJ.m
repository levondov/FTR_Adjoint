
verbose = false;
plotverbose = false;

fixedEpm = false;
useJvar = true;

% Initial conditions
X0 = [
    0.001306266362858
    -0.000681110843097
    0
    0.000082600114162
    0.000160506016655
    0
    0.121505236229312
    0.063375320588014
    0
    0
    0 ] * 1e-3;

% calculate E+ and E-
ex2 = 7.0e-6^2;
ey2 = 7.0e-6^2;
denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
X0(7) = Ep;
X0(8) = Em;

% create lattice
Xn = X0';

turns = 3;
numEvenCells = 5;

% setup moment object
mom = MomentSolverPeriodic(10e3, 0.0, X0);
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
dj0 = mom.CalcFoMGradientJ();

% gradient descent parameter
gamma = (f0/sum(df0.^2));
%%
% init arrays
gamma_h = gamma;
Xn_h = Xn;
f_h = f0;
fp_h = f0p;
df_h = df0;
dj_h = {dj0};

% adjust starting gamma
if useJvar
    Xn_h(end+1,:) = Xn - gamma_h(end)*df0' + gamma_h(end)*calcJcomponent(dj0, df0)'; % iterate
else
    Xn_h(end+1,:) = Xn - gamma_h(end)*df0';
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
%%
fprintf(['FoM: ',num2str(f_h(end)),'\n']);
while f_h(end) >= f0
    gamma_h(end+1) = gamma_h(end)/2.0;
    if useJvar
        Xn_h(end+1,:) = Xn - gamma_h(end)*df0' + gamma_h(end)*calcJcomponent(dj0, df0)'; % iterate
    else
        Xn_h(end+1,:) = Xn - gamma_h(end)*df0';
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
%%
if plotverbose
    figure('units','pixels','position',[100,100,1000,750]);
end
lastNpts = 100;
while 1
    ii=1;
    while f_h(end) < f_h(end-1)
        if plotverbose
            if (length(f_h) < lastNpts)
                subplot(3,3,1:3); cla();
                plot(log10(f_h));
                for i = 4:(4+length(an)-1)
                    subplot(3,3,i); cla();
                    plot(an_h(:,i-3));
                end
            else
                subplot(3,3,1:3); cla();
                plot(log10(f_h(end-100+1:end)));
                for i = 4:(4+length(an)-1)
                    subplot(3,3,i); cla();
                    plot(an_h(end-100+1:end,i-3));
                end
            end
            pause(0.01);
        end
        
        fprintf(['Iterating ',num2str(ii),'\n']);
        
        % iterate
        if useJvar
            Xn_h(end+1,:) = Xn_h(end,:) - gamma_h(end)*df_h(:,end)' + gamma_h(end)*calcJcomponent(dj_h{end}, df_h(:,end))';
        else
            Xn_h(end+1,:) = Xn_h(end,:) - gamma_h(end)*df_h(:,end)';
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
        
        % compute fom
        [f_h(end+1),fp_h(end+1,:)] = mom.GetFAndDF1();
        fprintf(['FoM: ',num2str(f_h(end)),'\n']);
        ii = ii + 1;
        
        %if (ii > 100) % if we have done more than 10 iterations without recalculating the gradient, increase gamma
        %    gamma_h(end+1) = gamma_h(end)*2.0;
        %end
    end
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    Xn_h(end+1,:) = Xn_h(end-1,:);
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
        dj = mom.CalcFoMGradientJ();
        dj_h{end+1} = dj;
    end
    
    if (ii == 2) % meaning no improving from recalculating gradient.
        % change gamma
        fprintf(['Updating Gamma \n']);
        f0n = f_h(end); Xnn = Xn_h(end,:);
        
        while f_h(end) >= f0n
            gamma_h(end+1) = gamma_h(end)/2.0;
            if useJvar
                Xn_h(end+1,:) = Xnn - gamma_h(end)*df_h(:,end)' + gamma_h(end)*calcJcomponent(dj_h{end}, df_h(:,end))'; % iterate
            else
                Xn_h(end+1,:) = Xn_h(end,:) - gamma_h(end)*df_h(:,end)';
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
    
    if f_h(end) < 1e-35
        break;
    end
    if length(f_h) > 100000
        break;
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










