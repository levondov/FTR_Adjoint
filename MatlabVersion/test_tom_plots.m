global k_perv k0 hardedge_flag norm_flag e1 e2
e1 = 1;
e2 = 1;
norm_flag = 0;
k0 = 10;
hardedge_flag = 1;

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

[z,y_nope,y_adj_nope,params,motion] = gd_cadj(a);
[O_nope,N_nope] = calcON2(z,y_nope,params,k_perv);
[z,ybase] = gd_y(a);
[fombase,~] = gd_F(a);

global k0 k_solv norm_flag
komegabase = k_solv(find(k_solv ~= 0,1));
NN = 20;
pert_s = linspace(0.50,1.50,NN);

%% perturb solenoid
fom_val1 = zeros(NN,1);
fom_int1 = zeros(NN,1);
x1 = zeros(NN,1);
a_0 = a;
[z,y,y_adj,params,motion] = gd_cadj(a_0);

for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_0(2) = a(2)*pert_s(i);
    pt1 = a2params(a_0);
    pt2 = a2params(a);
    komega = pt1(3)/rho;
    x1(i) = komega - komegabase;
    k_solv(1) = komega;
    ii = length(z);
    [FoM, FoMp] = get_F_and_dF(ybase,ii);
    fom_val1(i) = fombase - FoM;
    
    % change in FOM with integral and adjoint variables
    params = a2params(a_0);
    [O_pert,N_pert] = calcON2(z,y,params,k_perv);
    for j = 1:length(z)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end    
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);  
    ii = length(z);
    tmp2 =  k0^(-2) * ( y(ii,7) - 0.5*komega^2*y(ii,1) + k_perv ) * (-1) * y(ii,1) * komega *e1;
    tmp3 = ( y(ii,10) - komega*y(ii,1) ) * (-1) * y(ii,1) * e2;  
    fom_int1(i) = int_value + (tmp2 + tmp3) * (pt1(3)-pt2(3)) / rho;    
end

%%
figure; 
oset = 0;
hold on;
plot(-x1,fom_int1+oset); plot(x1,fom_val1); title('Solenoid strength'); 
legend('Integral','direct measurement','location','southeast');
xlabel('\Delta k_\Omega^2'); ylabel('\Delta FoM');
grid on;

%%
global k_perv k0 hardedge_flag e1 e2
k_perv = c2perv(0.0);
k0 = 10;
e1 = 0.0;
e2 = 0.0;
hardedge_flag = 1;

a = [
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1]';

[z,ybase] = gd_y(a);






