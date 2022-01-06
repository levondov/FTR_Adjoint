
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

%% Plot CoM across z 

CoM = mom.y(:,1).*mom.y(:,7) + mom.y(:,2).*mom.y(:,8) + mom.y(:,3).*mom.y(:,9) + mom.y(:,10).^2 - ...
    0.5 * (mom.y(:,4).^2 + mom.y(:,5).^2 + mom.y(:,6).^2);

xx = mom.y(:,1)+mom.y(:,2);
yy = mom.y(:,1)-mom.y(:,2);

figure;
subplot(1,2,1);
plot(mom.z,CoM); title('Constant of Motion');
subplot(1,2,2);
plot(mom.z, CoM ./ xx); title('Constant of Motion / <x^2>');

mom.PlotBeamEmit();

%% Calculate stuff after each iteration

[N,~] = size(Xn_h);

ex2 = zeros(N,1);
ey2 = zeros(N,1);
CoM = zeros(N,1);

idx = 1;
for i = 1:N
    fprintf([num2str(i),'/',num2str(N),'\n']);
    mom.initialMoments = Xn_h(i,:)';
    mom = mom.RunMoments(0);    
    ex2(i) = ( (mom.y(idx,7)+mom.y(idx,8)) .* (mom.y(idx,1)+mom.y(idx,2)) - 0.5*(mom.y(idx,4)+mom.y(idx,5)).^2 );
    ey2(i) = ( (mom.y(idx,7)-mom.y(idx,8)) .* (mom.y(idx,1)-mom.y(idx,2)) - 0.5*(mom.y(idx,4)-mom.y(idx,5)).^2 );
    CoM(i) = mom.y(idx,1).*mom.y(idx,7) + mom.y(idx,2).*mom.y(idx,8) + mom.y(idx,3).*mom.y(idx,9) + mom.y(idx,10).^2 - ...
        0.5 * (mom.y(idx,4).^2 + mom.y(idx,5).^2 + mom.y(idx,6).^2);    
end

%%
extrue = 7.6e-6^2;
figure; hold on;
plot(ex2); plot(ey2); legend('e_x^2','e_y^2');
yline(extrue,'linewidth',2); 
ylim([0,6]*1e-11);
figure; plot(CoM);


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

%%

figure; plot(Xn_h(:,[1,7]),'linewidth',2); legend('Q+','E+'); xlabel('Iterations'); ylabel('Moments'); grid on;
figure; plot(log10(fp_h),'linewidth',2); legend('|Q|^2','|P|^2','|E|^2','|L|^2'); xlabel('Iterations'); ylabel('log10(FoM)'); grid on;



