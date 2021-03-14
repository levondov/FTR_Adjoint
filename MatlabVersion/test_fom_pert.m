% Some test plots for perturbing magnets and measuring fom
params = [0.213,0.863,15e-4,...
    .00425+prof_offset,.0001*ql,-18.236*qs(1),...
    0.10655+prof_offset,0.0001*ql,21.3640*qs(2),...
    0.20895+prof_offset,0.0001*ql,-18.236*qs(3),...
    45*pi/180,45*pi/180,45*pi/180];
h=0.000001;
z_interval = [0.0,0.222]; % meters
z = z_interval(1):h:z_interval(2);
[y_nope] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
[O_nope,N_nope] = calcON2(z_adj,y_nope,params,k_perv);

%% perturb solenoid strength by +/- 20%
NN = 10;
pert_v = linspace(0.90,1.10,NN);
fom_val1 = zeros(NN,1);
fom_int1 = zeros(NN,1);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    params = [0.213,0.863,15e-4*pert_v(i),...
        .00425+prof_offset,.0001*ql,-18.236*qs(1),...
        0.10655+prof_offset,0.0001*ql,21.3640*qs(2),...
        0.20895+prof_offset,0.0001*ql,-18.236*qs(3),...
        45*pi/180,45*pi/180,45*pi/180];
    [y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    
    % left side of fom integration equation
    fom_val1(i) = (y(end,4:6)-y_nope(end,4:6))*init_cond1(4:6) - ...
        (y(end,1:3)-y_nope(end,1:3))*init_cond1(7:9) - ...
        (y(end,7:9)-y_nope(end,7:9))*init_cond1(1:3) - ...
        (y(end,10)-y_nope(end,10))*init_cond1(10);
    
    % right side of fom integration equation (to compare)
    [O_pert,N_pert] = calcON2(z_adj,y_nope,params,k_perv);
    for j = 1:length(z_adj)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end
    [int_value] = find_int(z_adj,y_adj,y_nope,O_pert,N_pert);
    fom_int1(i) = int_value;
end

%% perturb quad 1
pert_v = linspace(0.80,1.05,NN);
fom_val2 = zeros(NN,1);
fom_int2 = zeros(NN,1);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    params = [0.213,0.863,15e-4,...
        .00425+prof_offset,.0001*ql,-18.236*qs(1)*pert_v(i),...
        0.10655+prof_offset,0.0001*ql,21.3640*qs(2),...
        0.20895+prof_offset,0.0001*ql,-18.236*qs(3),...
        45*pi/180,45*pi/180,45*pi/180];
    [y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    
    % left side of fom integration equation
    fom_val2(i) = (y(end,4:6)-y_nope(end,4:6))*init_cond1(4:6) - ...
        (y(end,1:3)-y_nope(end,1:3))*init_cond1(7:9) - ...
        (y(end,7:9)-y_nope(end,7:9))*init_cond1(1:3) - ...
        (y(end,10)-y_nope(end,10))*init_cond1(10);
    
    % right side of fom integration equation (to compare)
    [O_pert,N_pert] = calcON2(z_adj,y_nope,params,k_perv);
    for j = 1:length(z_adj)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end
    [int_value] = find_int(z_adj,y_adj,y_nope,O_pert,N_pert);
    fom_int2(i) = int_value;
end

%% perturb quad 2
pert_v = linspace(0.90,1.10,NN);
fom_val3 = zeros(NN,1);
fom_int3 = zeros(NN,1);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    params = [0.213,0.863,15e-4,...
        .00425+prof_offset,.0001*ql,-18.236*qs(1),...
        0.10655+prof_offset,0.0001*ql,21.3640*qs(2)*pert_v(i),...
        0.20895+prof_offset,0.0001*ql,-18.236*qs(3),...
        45*pi/180,45*pi/180,45*pi/180];
    [y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    
        % left side of fom integration equation
    fom_val3(i) = (y(end,4:6)-y_nope(end,4:6))*init_cond1(4:6) - ...
        (y(end,1:3)-y_nope(end,1:3))*init_cond1(7:9) - ...
        (y(end,7:9)-y_nope(end,7:9))*init_cond1(1:3) - ...
        (y(end,10)-y_nope(end,10))*init_cond1(10);
    
    % right side of fom integration equation (to compare)
    [O_pert,N_pert] = calcON2(z_adj,y_nope,params,k_perv);
    for j = 1:length(z_adj)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end
    [int_value] = find_int(z_adj,y_adj,y_nope,O_pert,N_pert);
    fom_int3(i) = int_value;
end


%% perturb quad 3
fom_val4 = zeros(NN,1);
fom_int4 = zeros(NN,1);
for i = 1:NN
    fprintf([num2str(i),'/',num2str(NN),'\n']);
    params = [0.213,0.863,15e-4,...
        .00425+prof_offset,.0001*ql,-18.236*qs(1),...
        0.10655+prof_offset,0.0001*ql,21.3640*qs(2),...
        0.20895+prof_offset,0.0001*ql,-18.236*qs(3)*pert_v(i),...
        45*pi/180,45*pi/180,45*pi/180];
    [y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    
            % left side of fom integration equation
    fom_val4(i) = (y(end,4:6)-y_nope(end,4:6))*init_cond1(4:6) - ...
        (y(end,1:3)-y_nope(end,1:3))*init_cond1(7:9) - ...
        (y(end,7:9)-y_nope(end,7:9))*init_cond1(1:3) - ...
        (y(end,10)-y_nope(end,10))*init_cond1(10);
    
    % right side of fom integration equation (to compare)
    [O_pert,N_pert] = calcON2(z_adj,y_nope,params,k_perv);
    for j = 1:length(z_adj)
        O_pert{j} = O_pert{j} - O_nope{j};
        N_pert{j} = N_pert{j} - N_nope{j};
    end
    [int_value] = find_int(z_adj,y_adj,y_nope,O_pert,N_pert);
    fom_int4(i) = int_value;
end


%%
figure; 

subplot(2,2,1); hold on;
plot(pert_v,fom_int1); plot(pert_v,fom_val1); title('Solenoid strength'); legend('Integral','direct measurement');

subplot(2,2,2); hold on;
plot(pert_v,fom_int2); plot(pert_v,fom_val2); title('quad 1 strength');

subplot(2,2,3); hold on;
plot(pert_v,fom_int3); plot(pert_v,fom_val3); title('quad 2 strength');

subplot(2,2,4); hold on;
plot(pert_v,fom_int4); plot(pert_v,fom_val4); title('quad 3 strength');

stop
%%
figure;

subplot(4,3,1);
plot(pert_v,fom_val1(:,1:3),'.-'); legend({'Q+','Q-','Qx'});

subplot(4,3,2);
plot(pert_v,fom_val1(:,4:6),'.-'); legend({'P+','P-','Px'});

subplot(4,3,3);
plot(pert_v,fom_val1(:,7:9),'.-'); legend({'E+','E-','Ex'});

%%%%%
subplot(4,3,4);
plot(pert_v,fom_val2(:,1:3),'.-'); legend({'Q+','Q-','Qx'});

subplot(4,3,5);
plot(pert_v,fom_val2(:,4:6),'.-'); legend({'P+','P-','Px'});

subplot(4,3,6);
plot(pert_v,fom_val2(:,7:9),'.-'); legend({'E+','E-','Ex'});

%%%%%
subplot(4,3,7);
plot(pert_v,fom_val3(:,1:3),'.-'); legend({'Q+','Q-','Qx'});

subplot(4,3,8);
plot(pert_v,fom_val3(:,4:6),'.-'); legend({'P+','P-','Px'});

subplot(4,3,9);
plot(pert_v,fom_val3(:,7:9),'.-'); legend({'E+','E-','Ex'});

%%%%%
subplot(4,3,10);
plot(pert_v,fom_val4(:,1:3),'.-'); legend({'Q+','Q-','Qx'});

subplot(4,3,11);
plot(pert_v,fom_val4(:,4:6),'.-'); legend({'P+','P-','Px'});

subplot(4,3,12);
plot(pert_v,fom_val4(:,7:9),'.-'); legend({'E+','E-','Ex'});




