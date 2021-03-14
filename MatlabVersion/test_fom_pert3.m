global k_perv k0 hardedge_flag norm_flag
norm_flag = 1;
k0 = 10;
hardedge_flag = 1;

a = [
    1.1646
    1.0182
    0.9993
    1.0100
    1.0157
    1.0381
    1.0102
    1.0380
    1.0251
    1.0278
    1.0155]';

[z,y_nope,y_adj_nope,params,motion] = gd_cadj(a);
[O_nope,N_nope] = calcON2(z,y_nope,params,k_perv);
fom_base = gd_F(a);

global k0 k_solvv norm_flag
norm_flag = 1
NN = 20;
pert_s = linspace(0.50,1.50,NN);
pert_v = linspace(0.75,1.25,NN);
pert_v2 = linspace(0.90,1.1,NN);
pert_v3 = linspace(0.75,1.25,NN);

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
    x1(i) = pt1(3)^2-pt2(3)^2;
    %x1(i) = a_0(2)-a(2);
    
    % direct FOM calculation    
    fom_val1(i) = fom_base - gd_F(a_0);
    
    % change in FOM with integral and adjoint variables
    params = a2params(a_0);
    [O_pert,N_pert] = calcON2(z,y,params,k_perv);
    for j = 1:length(z)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end    
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);
    fom_int1(i) = int_value;
end

%% perturb quad 1
fom_val2 = zeros(NN,1);
fom_int2 = zeros(NN,1);
x2 = zeros(NN,1);
a_1 = a;
[z,y,y_adj,params,motion] = gd_cadj(a_1);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_1(4) = a(4)*pert_v(i);
    pt1 = a2params(a_1);
    pt2 = a2params(a);    
    x2(i) = pt1(6) - pt2(6);
    %x2(i) = a_1(4)-a(4);
    
    % direct FOM calculation    
    fom_val2(i) = fom_base - gd_F(a_1);
    
    % change in FOM with integral and adjoint variables
    params = a2params(a_1);
    [O_pert,N_pert] = calcON2(z,y,params,k_perv);
    for j = 1:length(z)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end    
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);
    fom_int2(i) = int_value;
end

%% perturb quad 2
fom_val3 = zeros(NN,1);
fom_int3 = zeros(NN,1);
x3 = zeros(NN,1);
a_2 = a;
[z,y,y_adj,params,motion] = gd_cadj(a_2);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_2(6) = a(6)*pert_v2(i);
    pt1 = a2params(a_2);
    pt2 = a2params(a);    
    x3(i) = pt1(9)-pt2(9);
    %x3(i) = a_2(6)-a(6);
    
    % direct FOM calculation    
    fom_val3(i) = fom_base - gd_F(a_2);
    
    % change in FOM with integral and adjoint variables
    params = a2params(a_2);
    [O_pert,N_pert] = calcON2(z,y,params,k_perv);
    for j = 1:length(z)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end    
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);
    fom_int3(i) = int_value;
end

%% perturb quad 3
fom_val4 = zeros(NN,1);
fom_int4 = zeros(NN,1);
x4 = zeros(NN,1);
a_3 = a;
[z,y,y_adj,params,motion] = gd_cadj(a_3);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    a_3(8) = a(8)*pert_v3(i);
    pt1 = a2params(a_3);
    pt2 = a2params(a);    
    x4(i) = pt1(12)-pt2(12);
    %x4(i) = a_3(8)-a(8);
    
    % direct FOM calculation    
    fom_val4(i) = fom_base - gd_F(a_3);
    
    % change in FOM with integral and adjoint variables
    params = a2params(a_3);
    [O_pert,N_pert] = calcON2(z,y,params,k_perv);
    for j = 1:length(z)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end    
    [int_value] = find_int(z,y_adj,y,O_pert,N_pert);
    fom_int4(i) = int_value;
end

save('data_paper_1213/gradcomp/pert_st2.mat','fom_int1','fom_int2','fom_int3','fom_int4',...
    'fom_val1','fom_val2','fom_val3','fom_val4',...
    'pert_v','pert_s','x1','x2','x3','x4')


%%
if 0
   fom_val1 = abs(fom_val1); 
   fom_val2 = abs(fom_val2); 
   fom_val3 = abs(fom_val3); 
   fom_val4 = abs(fom_val4);   
   fom_int1 = abs(fom_int1);
   fom_int2 = abs(fom_int2);
   fom_int3 = abs(fom_int3);
   fom_int4 = abs(fom_int4);   
end
figure; 
oset = 0;
subplot(2,2,1); hold on;
plot(-x1,fom_int1+oset); plot(x1,fom_val1); title('Solenoid strength'); 
legend('Integral','direct measurement','location','southeast');
grid on;

subplot(2,2,2); hold on;
plot(-x2,fom_int2+oset); plot(x2,fom_val2); title('quad 1 strength');
grid on;

subplot(2,2,3); hold on;
plot(-x3,fom_int3+oset); plot(x3,fom_val3); title('quad 2 strength');
grid on;

subplot(2,2,4); hold on;
plot(-x4,fom_int4+oset); plot(x4,fom_val4); title('quad 3 strength');
grid on;

%figure; hold on;
%plot(fom_int1+oset);
%plot(fom_int2+oset);
%plot(fom_int3+oset);
%plot(fom_int4+oset);