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
               
momX = MomentSolverPeriodicChaser(10e3, 0, X0); momX.h = 100; % default is h=10000 which is too much for now

% Create our lattice of 16 periods 
% here we have 32 quads each with a adjustable strength dB, thus our a vector is length 32
period = 1;
qlength = 0.02;
qstart = [0.0, 0.03, 0.07];
db = [0.5, -0.5, 0.5];
dbfull = repmat(db,1,period);
a = repmat(ones(1,length(db)),1,period);
qend = [qstart(1)+qlength/2.0, qstart(2)+qlength, qstart(3)+qlength/2.0];
qrot = [0.1, 0.1, 0.1];
zstart = 0.0;
zend = 0.08;

% setup lattice
momX = momX.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

% run moment + adjoint equations
momX = momX.RunMoments();
momX = momX.RunMomentsAdjoint();

O_nope = momX.Omat;
N_nope = momX.Nmat;
FoM_base = momX.GetW;
lattice_base = momX.lattice;

NN = 40;
pert_db = linspace(0.95,1.05,NN);

% perturb db 1
fom_val1 = zeros(NN,1);
fom_int1 = zeros(NN,1);
x1 = zeros(NN,1);
a_0 = a;
idx = 1;
momXbase = momX;
for i = 1:NN   
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    % perturb a 
    a_0(idx) = a(idx)*pert_db(i);    
    x1(i) = dbfull(idx)*a_0(idx) - dbfull(idx);
    
    % direct FOM calculation
    momX.lattice = lattice_base;
    momX = momX.UpdateLatticeProfile(a_0);
    momX = momX.RunMoments();
    fom_val1(i) = FoM_base - momX.GetW;
    
    % change in FOM with integral and adjoint variables
    [O_pert, N_pert, ~] = momXbase.CalcONandACTMatrix(momX.lattice);
    for j = 1:length(momXbase.z)
        O_pert{j} = (O_pert{j} - O_nope{j});
        N_pert{j} = (N_pert{j} - N_nope{j});
    end  
    [int_value] = momXbase.AdjointIntegral(O_pert, N_pert);
    fom_int1(i) = int_value;    
end

% plot

figure; 
subplot(1,3,1);
oset = 0;
 hold on;
plot(-x1,fom_int1+oset,'.-'); plot(x1,fom_val1,'.-'); title('Quad strength dFE/dX'); 
legend('Integral','direct measurement','location','southeast');
grid on;

subplot(1,3,2);
plot(-x1,fom_int1,'.-')
grid on;
legend('Integral alone');

subplot(1,3,3); hold on;
plot(0,0);
plot(x1,fom_val1,'.-')
grid on;
legend('','direct alone');

%% perturb FE
NN = 40;
pert_db = linspace(0.7,1.6,NN);

% setup lattice
momX = MomentSolverPeriodicChaser(10e3, 0, X0); momX.h = 100; % default is h=10000 which is too much for now
momX = momX.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

% run moment + adjoint equations
momX.M = 1;
momX.QTE = 1e-7;
momX = momX.RunMoments();
momX = momX.RunMomentsAdjointFE();

O_nope = momX.Omat;
N_nope = momX.Nmat;
FoM_base = momX.GetFE;
lattice_base = momX.lattice;

fom_val1 = zeros(NN,1);
fom_int1 = zeros(NN,1);
x1 = zeros(NN,1);
a_0 = a;
idx = 1;

momXbase = momX;
for i = 1:NN  
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    % perturb a 
    a_0(idx) = a(idx)*pert_db(i);    
    x1(i) = dbfull(idx)*a_0(idx) - dbfull(idx);
    
    % direct FOM calculation
    momX.lattice = lattice_base;
    momX = momX.UpdateLatticeProfile(a_0);
    momX = momX.RunMoments();
    fom_val1(i) = (FoM_base - momX.GetFE);
    
    % change in FOM with integral and adjoint variables
    [O_pert, N_pert, ~] = momXbase.CalcONandACTMatrix(momX.lattice);
    for j = 1:length(momXbase.z)
        O_pert{j} = (O_pert{j} - O_nope{j});
        N_pert{j} = (N_pert{j} - N_nope{j});
    end  
    [int_value] = momXbase.AdjointIntegral(O_pert, N_pert);
    fom_int1(i) = -int_value;    
end

figure; 
subplot(1,3,1);
oset = 0;
 hold on;
plot(-x1,fom_int1+oset,'.-'); plot(x1,fom_val1,'.-'); title('Quad strength dFE/da'); 
legend('Integral','direct measurement','location','southeast');
grid on;

subplot(1,3,2);
plot(-x1,fom_int1,'.-')
grid on;
legend('Integral alone');

subplot(1,3,3); hold on;
plot(0,0);
plot(x1,fom_val1,'.-')
grid on;
legend('','direct alone');

%% Perturb X variables
% setup lattice
momXbase = MomentSolverPeriodicChaser(10e3, 0, X0); momXbase.h = 100; % default is h=10000 which is too much for now
momXbase = momXbase.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

% run moment + adjoint equations
momXbase = momXbase.RunMoments();
momXbase = momXbase.RunMomentsAdjoint();
df0 = momXbase.GetWGradientX();

FoM_base = momXbase.GetW;
X0_nope = X0;

NN = 40;
pert_db = linspace(0.95,1.05,NN);

% perturb db 1
fom_val1 = zeros(NN,1);
fom_int1 = zeros(NN,1);
x1 = zeros(NN,1);
Xp = X0;
idx = 1;
momX = momXbase;
baseVal = [momXbase.y(1,idx), momXbase.y(end,idx)];
pert_db(21)=1.0;
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    % perturb X
    Xp(idx) = X0(idx) * pert_db(i);
    x1(i) = pert_db(i);  
    
    % direct FOM calculation
    momX.initialMoments = Xp;
    momX = momX.RunMoments();
    fom_val1(i) = (FoM_base - momX.GetW);
end
fom_int1 = momX.k0 * X0(idx) * momXbase.GetWGradientX * (1 - pert_db);

% plot

figure; 
subplot(1,3,1);
oset = 0;
 hold on;
plot(x1,fom_int1(idx,:)+oset,'.-'); plot(x1,fom_val1,'.-'); title('Quad strength dQ+/dX'); 
legend('adjoint','direct measurement','location','southeast');
grid on;

subplot(1,3,2);
plot(x1,fom_int1(idx,:),'.-')
grid on;
legend('adjoint alone');

subplot(1,3,3); hold on;
plot(1,0);
plot(x1,fom_val1,'.-')
grid on;
legend('','direct alone');