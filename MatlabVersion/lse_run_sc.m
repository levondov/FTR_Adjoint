%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization , only run once
%
% initial seed and distributions
rng(1234);
[init_cond] = lse_init;

% setup flags for the solver
global hardedge_flag k_perv
k_perv = c2perv(1e-3);
hardedge_flag = 1;
% some extra parameters we want to calculate during the simulations
global extra_params_flag extra_params_flag_firstcall
extra_params_flag = 1;
extra_params_flag_firstcall = 1;

% Initial magnet and drift positions
% [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
ql = 1;
qs = [1,1,1];
prof_offset = 0;
params = [0.213,0.863,15e-4,.00425+prof_offset,.0001*ql,-18.236*qs(1),0.10655+prof_offset,0.0001*ql,...
    21.3640*qs(2),0.20895+prof_offset,0.0001*ql,-18.236*qs(3),45*pi/180,45*pi/180,45*pi/180];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% run integration to solve diff eqn.
h=0.000001;
z_interval = [0.0,0.212]; % meters
z = z_interval(1):h:z_interval(2);
tic
[y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
toc

% Constant of Motion %%
L = [y(:,10)];
EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
motion = EQ + (1/2)*L.^2 - (1/2)*PP;

dphi = diff(y(:,11))./diff(z');
dphi(end+1) = dphi(end);
%y = lar2cart(y,dphi);
lseplot(z,y,motion,'');

%%
perv = [c2perv(linspace(0,0.1e-3,50))];
extra_params_flag = 0; extra_params_flag_firstcall = 1;
% run integration to solve diff eqn.
h=0.000001;
z_interval = [0.0,0.212]; % meters
z = z_interval(1):h:z_interval(2);

yy = cell(length(perv),1);
for i = 1:length(perv)
    fprintf(['Run ',num2str(i),'/',num2str(length(perv)),'\n']);
    k_perv = perv(i);
    [y] = ode3(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    yy{i} = y(1:100:end,:);
end
z=z(1:100:end);

%%
idx_r = 25; %lenght(perv);
pp_points=1;
clrs = jet(idx_r);
titles = {'Q+','Q-','Qx','P+','P-','Px','E+','E-','Ex'};

figure;
for ii = 1:9
    subplot(3,3,ii); hold on;
    for i = 1:idx_r
        plot(z,yy{i}(1:pp_points:end,ii),'color',clrs(i,:));
    end
    title(titles{ii});
    %set(gca, 'YScale', 'log');
end

hp4 = get(subplot(3,3,9),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025  hp4(2)+hp4(3)*3.2])
colormap(jet); caxis([perv(1),perv(idx_r)])

%% Just plot final point.
%crts = [linspace(0,0.1e-3,50)];
yy_finalpoints = zeros(length(yy),10);

for i = 1:length(yy)
   yy_finalpoints(i,:) = yy{i}(end,1:10);    
end

idx_r = 50;
splot_idx = [1,1,1,2,2,2,3,3,3,4];
figure;
for i = 1:10
    subplot(2,2,splot_idx(i)); hold on;
    val = 2*abs(yy_finalpoints(1,i) - yy_finalpoints(1:idx_r,i))./(yy_finalpoints(1,i)+yy_finalpoints(1:idx_r,i));
    %val = 2*abs(yy_finalpoints(1,i) - yy_finalpoints(1:idx_r,i))./(yy_finalpoints(1,i));
    plot(crts(1:idx_r)*1e3,val,'.-');
    xlabel('Current (mA)');
end

subplot(2,2,1);
legend('Q+','Q-','Qx');
subplot(2,2,2);
legend('P+','P-','Px');
subplot(2,2,3);
legend('E+','E-','Ex');
subplot(2,2,4);
legend('L');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'fractional percentage change from 0 current solution at solenoid entrance', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize',12)





