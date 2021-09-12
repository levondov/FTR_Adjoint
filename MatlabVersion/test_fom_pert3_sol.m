global k_perv k0 hardedge_flag norm_flag e1 e2
e1 = 1;
e2 = 1;
norm_flag = 0;
k0 = 10;
hardedge_flag = 1;
k_perv = c2perv(0.0);

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);

a = [
   1.159092901897256
  -0.979394946584468
   0.999359301160531
   1.010646218457928
   1.013283550954570
   1.041899091927412
   1.006294655428367
   1.050000353512966
   1.027780261173520
   1.036071870631294
   1.022099289148964]';

global k0 k_solvv norm_flag komegaE
NN = 20;
pert_s = linspace(0.95,1.05,NN);

%% perturb solenoid
fom_val1 = zeros(NN,1);
fom_val2 = zeros(NN,1);
fom_val3 = zeros(NN,1);
fom_int1 = zeros(NN,1);
fom_int2 = zeros(NN,1);
x1 = zeros(NN,1);
a_0 = a;
pt2 = a2params(a);

% reference case
komegaE = pt2(3) / rho;
fom_base = gd_F(a);

% calculate implicit change
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_0(2) = a(2)*pert_s(i);
    pt1 = a2params(a_0);
    x1(i) = (pt1(3)^2-pt2(3)^2) / rho^2;
    
    komegaE = pt2(3) / rho;
    
    % direct FOM calculation    
    fom_val1(i) = (gd_F(a_0) - fom_base);
end

% calculate explicit change
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_0(2) = a(2)*pert_s(i);
    pt1 = a2params(a_0);
    x1(i) = (pt1(3)^2-pt2(3)^2) / rho^2;
    
    komegaE = pt1(3) / rho;
    
    % direct FOM calculation    
    fom_val2(i) = (gd_F(a) - fom_base);
end

% calculate total change
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_0(2) = a(2)*pert_s(i);
    pt1 = a2params(a_0);
    x1(i) = (pt1(3)^2-pt2(3)^2) / rho^2;
    
    komegaE = pt1(3) / rho;
    
    % direct FOM calculation    
    fom_val3(i) = (gd_F(a_0) - fom_base);
end

%%
% reference adjoint
komegaE = pt2(3) / rho;
[z,y,y_adj,params,motion] = gd_cadj(a);
[O_nope,N_nope] = calcON2(z,y,params,k_perv);

% calculate change of O,N matrix
params2 = params;
params2(3) = params(3) + params(3)*0.01;
k1 = params(3) / rho; k2 = params2(3) / rho;
[O_pert,N_pert] = calcON2(z,y,params2,k_perv);
for j = 1:length(z)
    O_pert{j} = (O_pert{j} - O_nope{j}) / (k2^2 - k1^2);
    N_pert{j} = (N_pert{j} - N_nope{j}) / (k2^2 - k1^2);
end

[z,y,m1,m2,params] = gd_y(a);

% calculate adjoint calculations
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_0(2) = a(2)*pert_s(i);
    pt1 = a2params(a_0);
    x1(i) = (pt1(3)^2-pt2(3)^2) / rho^2;

    % change in FOM with integral and adjoint variables
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);
    
    % implcit change
    fom_int1(i) = x1(i) * int_value;
    
    % explicit change
    ii = length(z);
    %tmp2 =  k0^(-2) * ( y(ii,7) - 0.5*komega^2*y(ii,1) + k_perv ) * (-1) * y(ii,1) * komega *e1;
    f5_tmp = ( y(ii,7) + 0.5*komegaE^2*y(ii,1) - komegaE*y(ii,10) );
    tmp3 = k0^(-2) * f5_tmp * e2 * ( 0.5*y(ii,1) - 0.5*y(ii,10)*(1/komegaE) );  
    fom_int2(i) = x1(i) * (tmp3);
end

%%
figure; 
hold on;
plot(x1,fom_val1,'linewidth',2,'Color',[0, 0.4470, 0.7410]);
plot(x1,fom_val2,'linewidth',2,'Color',	[0.8500, 0.3250, 0.0980]);
plot(x1,fom_val1+fom_val2,'linewidth',2, 'Color', [0.9290, 0.6940, 0.1250]);
plot(x1,fom_val3,'linewidth',2, 'Color', [0.4940, 0.1840, 0.5560]);
plot(x1,fom_val1+fom_val2-fom_val3,'linewidth',2,'Color','g')

plot(x1,fom_int1,'linewidth',2,'linestyle','--','Color',	[0, 0.4470, 0.7410])
plot(x1,fom_int2,'linewidth',2,'linestyle','--','Color',	[0.8500, 0.3250, 0.0980])

legend('\Delta F_1','\Delta F_2','\Delta F_1 + \Delta F_2','\Delta F_3','\Delta F_1 + \Delta F_2 - \Delta F_3','\Delta F_1^{(y)}','\Delta F_2^{(y)}','location','southeast');

title('Solenoid strength'); 
xlabel('\Delta k_\Omega^2'); ylabel('\Delta FoM');
grid on;

%%










