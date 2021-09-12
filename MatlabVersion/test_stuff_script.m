figure; subplot(2,2,1); hold on;
plot(z,y(:,1),'k');
plot(z,y(:,2),'k--');
plot(z,y(:,3),'k.-');
plot(z2,y2(:,1),'r')
plot(z2,y2(:,2),'r--')
plot(z2,y2(:,3),'r.-')
legend('Q+','Q-','Qx')

subplot(2,2,2); hold on;
plot(z,y(:,4),'k'); 
plot(z,y(:,5),'k--');
plot(z,y(:,6),'k.-'); 
plot(z2,y2(:,4),'r');
plot(z2,y2(:,5),'r--');
plot(z2,y2(:,6),'r.-');
legend('P+','P-','Px')

subplot(2,2,3); hold on;
plot(z,y(:,7),'k'); 
plot(z,y(:,8),'k--'); 
plot(z,y(:,9),'k.-'); 
plot(z2,y2(:,7),'r');
plot(z2,y2(:,8),'r--');
plot(z2,y2(:,9),'r.-');
legend('E+','E-','Ex')

subplot(2,2,4); hold on;
plot(z,y(:,10),'k--');
plot(z,y2(:,10),'r--');
legend('L');

figure; hold on;
plot(z,abs(y(:,1)-y2(:,1)))
plot(z,abs(y(:,2)-y2(:,2)))
plot(z,abs(y(:,3)-y2(:,3)))
plot(z,abs(y(:,4)-y2(:,4)))
plot(z,abs(y(:,5)-y2(:,5)))
plot(z,abs(y(:,6)-y2(:,6)))
plot(z,abs(y(:,7)-y2(:,7)))
plot(z,abs(y(:,8)-y2(:,8)))
plot(z,abs(y(:,9)-y2(:,9)))
plot(z,abs(y(:,10)-y2(:,10)))
set(gca, 'YScale', 'log')

figure;
subplot(2,2,1); hold on;
for i = 1:3
    yerr = 2*(y(:,i)-y2(:,i))./(abs(y(:,i))+abs(y2(:,i)));
    plot(z,yerr,'linewidth',2);
end
legend({'Q+','Q-','Qx'}); xlabel('Z position'); ylabel('Fractional Percent'); title('Difference between 0 and 1mA beam');
subplot(2,2,2); hold on;
for i = 4:6
    yerr = 2*(y(:,i)-y2(:,i))./(abs(y(:,i))+abs(y2(:,i)));
    plot(z,yerr,'linewidth',2);
end
legend({'P+','P-','Px'});
subplot(2,2,3); hold on;
for i = 7:9
    yerr = 2*(y(:,i)-y2(:,i))./(abs(y(:,i))+abs(y2(:,i)));
    plot(z,yerr,'linewidth',2);
end
legend({'E+','E-','Ex'});
subplot(2,2,4); hold on;
i=10;
yerr = 2*(y(:,i)-y2(:,i))./(abs(y(:,i))+abs(y2(:,i)));
plot(z,yerr,'linewidth',2);
legend({'L'});


%%
figure; subplot(2,2,1); hold on;
plot(z,y(:,1),'k'); 
plot(z,y(:,2),'k--');
plot(z,y(:,3),'k.-'); 
plot(z2,flip(y2(:,1)),'r');
plot(z2,flip(y2(:,2)),'r--');
plot(z2,flip(y2(:,3)),'r.-');
legend('Q+','Q-','Qx')
subplot(2,2,2); hold on;
plot(z,y(:,4),'k');
plot(z,y(:,5),'k--');
plot(z,y(:,6),'k.-');
plot(z2,flip(y2(:,4)),'r')
plot(z2,flip(y2(:,5)),'r--')
 plot(z2,flip(y2(:,6)),'r.-')
legend('P+','P-','Px')
subplot(2,2,3); hold on;
plot(z,y(:,7),'k');
plot(z,y(:,8),'k--');
plot(z,y(:,9),'k.-');
legend('E+','E-','Ex')
plot(z2,flip(y2(:,7)),'r')
plot(z2,flip(y2(:,8)),'r--')
plot(z2,flip(y2(:,9)),'r.-')
subplot(2,2,4); hold on;
plot(z,y(:,10),'k--'); hold on; plot(z2,flip(y2(:,10)),'r--')
legend('L');

figure; hold on;
plot(z,abs(y(:,1)-flip(y2(:,1))))
plot(z,abs(y(:,2)-flip(y2(:,2))))
plot(z,abs(y(:,3)-flip(y2(:,3))))
plot(z,abs(y(:,4)-flip(y2(:,4))))
plot(z,abs(y(:,5)-flip(y2(:,5))))
plot(z,abs(y(:,6)-flip(y2(:,6))))
plot(z,abs(y(:,7)-flip(y2(:,7))))
plot(z,abs(y(:,8)-flip(y2(:,8))))
plot(z,abs(y(:,9)-flip(y2(:,9))))
plot(z,abs(y(:,10)-flip(y2(:,10))))
set(gca, 'YScale', 'log')




%%
params(9) = -18;
[O1,N1] = calcON2(z_adj,y,params,k_perv);
params(9) = -20;
[O2,N2] = calcON2(z_adj,y,params,k_perv);

tmp=zeros(length(z_adj),9);
for j = 1:length(z_adj)
    O2{j} = O2{j} - O1{j};
    N2{j} = N2{j} - N1{j};
    tmp(j,1)=O2{j}(1,1);
    tmp(j,2)=O2{j}(1,2);
    tmp(j,3)=O2{j}(1,3);
    tmp(j,4)=O2{j}(2,1);
    tmp(j,5)=O2{j}(2,2);
    tmp(j,6)=O2{j}(2,3);
    tmp(j,7)=O2{j}(3,1);
    tmp(j,8)=O2{j}(3,2);
    tmp(j,9)=O2{j}(3,3);    
end
figure; 
subplot(3,3,1); plot(tmp(:,1));
subplot(3,3,2); plot(tmp(:,2));
subplot(3,3,3); plot(tmp(:,3));
subplot(3,3,4); plot(tmp(:,4));
subplot(3,3,5); plot(tmp(:,5));
subplot(3,3,6); plot(tmp(:,6));
subplot(3,3,7); plot(tmp(:,7));
subplot(3,3,8); plot(tmp(:,8));
subplot(3,3,9); plot(tmp(:,9));


%%%%%%
%%

load('opt_2_0mA_norm.mat');
ff=zeros(length(f_h_new),1);
jj=1;
ff(jj) = f_h_new(1);

for i = 1:length(f_h_new)
    if f_h_new(i) < ff(jj)
       ff(jj+1) = f_h_new(i);
       jj=jj+1;
    end
end
ff1=ff(1:find(ff==0,1)-1);

load('opt_3_1mA_norm.mat');
ff=zeros(length(f_h_new),1);
jj=1;
ff(jj) = f_h_new(1);

for i = 1:length(f_h_new)
    if f_h_new(i) < ff(jj)
       ff(jj+1) = f_h_new(i);
       jj=jj+1;
    end
end
ff2=ff(1:find(ff==0,1)-1);
%%

load('0f.mat');
ff=zeros(length(f_h_new),1);
jj=1;
ff(jj) = f_h_new(1);

for i = 1:length(f_h_new)
    if f_h_new(i) < ff(jj)
       ff(jj+1) = f_h_new(i);
       jj=jj+1;
    end
end
ff3=ff(1:find(ff==0,1)-1);

load('0v.mat');
ff=zeros(length(f_h_new),1);
jj=1;
ff(jj) = f_h_new(1);

for i = 1:length(f_h_new)
    if f_h_new(i) < ff(jj)
       ff(jj+1) = f_h_new(i);
       jj=jj+1;
    end
end
ff4=ff(1:find(ff==0,1)-1);

%%

figure; 

subplot(1,2,1); hold on;
%plot(ff1(10:end),'linewidth',2);
plot(ff2(10:end),'linewidth',2);
plot(ff3(10:end),'linewidth',2);
plot(ff4(10:end),'linewidth',2);
grid on; xlabel('Iteration'); ylabel('FoM / k0^2 Q_+^2(0)');

subplot(1,2,2); hold on;
%plot(log(ff1(10:end)),'linewidth',2);
plot(log10(ff2(10:end)),'linewidth',2);
plot(log10(ff3(10:end)),'linewidth',2);
plot(log10(ff4(10:end)),'linewidth',2);
grid on; xlabel('Iteration'); ylabel('log(FoM / k0^2 Q_+^2(0))');


figure; hold on;
plot(log10(abs(diff(ff2(10:end)))/ff2(10)),'linewidth',2);
plot(log10(abs(diff(ff3(10:end)))/ff3(10)),'linewidth',2);
plot(log10(abs(diff(ff4(10:end)))/ff4(10)),'linewidth',2);

xlabel('Iteration');
ylabel('log((FoM_{n+1} - FoM_n) / FoM_0)');



%%
yy=y_old;
idx=24000;
x2 = y(idx:end,1)+y(idx:end,2);
y2 = y(idx:end,1)-y(idx:end,2);
figure; hold on;
plot(z(idx:end),x2);
plot(z(idx:end),y2);


%%

[z,y,m1,m2] = gd_y(an_h(end,:));
lseplot(z,y,m1,'')

figure; 
plot(log10(fp_h),'linewidth',3); legend('FoM 1','FoM 2','FoM 3','FoM 4','FoM 5')
title(['k0 = ',num2str(k0),' | e1 = ',num2str(e1),' | e2 = ',num2str(e2),' | \Lambda = ',num2str(k_perv)]);
ylabel('Log10(FoM)'); xlabel('Iterations');

global zv k_solv
komega = max(k_solv);
figure;
subplot(2,1,1);
[z,y,m1,m2] = gd_y(an_h(1,:));
hold on;
plot(z(1:24491), -y((1:24491),1),'linewidth',2);
plot(z(24491:end),-y((24491:end),1)*komega,'linewidth',2);
plot(z,y(:,10),'linewidth',2);
plot(zv,k_solv*1e-6-3e-5,'linewidth',2);
grid on;
legend('-Q_+','-Q_+ k_\Omega','L','1e-6k_\Omega - 3e-5')
title('FoM 5 components pre optimization');
xlabel('Z position');
subplot(2,1,2);
[z,y,m1,m2] = gd_y(an_h(end,:));
hold on;
plot(z(1:24491), -y((1:24491),1),'linewidth',2);
plot(z(24491:end),-y((24491:end),1)*komega,'linewidth',2);
plot(z,y(:,10),'linewidth',2);
plot(zv,k_solv*1e-6-3e-5,'linewidth',2);
grid on;
legend('-Q_+','-Q_+ k_\Omega','L','1e-6k_\Omega - 3e-5')
title('FoM 5 components post optimization');
xlabel('Z position');


%%

figure; hold on;
plot(log10(fp_h(:,1:4)),'linewidth',3); 
plot(log10(fp_h(:,5)*1e10),'linewidth',3);
legend('FoM 1','FoM 2','FoM 3','FoM 4','1e10*FoM 5')
title(['k0 = ',num2str(k0),' | \Lambda = ',num2str(k_perv)]);
ylabel('Log10(FoM)'); xlabel('Iterations');

%%

kval = zeros(length(an_h(:,1)),1);
for i = 1:length(kval)
   tmp = a2params(an_h(i,:));
   kval(i) = tmp(3) / rho;
end

figure; hold on;

yyaxis left;
plot(10:length(kval),kval(10:end))
%ylim([-1,-0.975])

yyaxis right;
plot(10:length(kval),log10(fp_h(10:end,[1,2,3])));
%ylim([-17,-8]);

legend('k_\Omega','FoM1','FoM2','FoM3');
%legend('k_\Omega','k_0^{-2} (E_+ - 0.5 k_\Omega^2 Q_+ + \Lambda)^2')
%legend('a_3','a_4','a_5','a_6','a_7','a_8','a_9','a_10','a_11','k_0^{-2} (E_+ + 0.5 k_\Omega^2 Q_+ - k_\Omega L)^2','Location','southeast')
%legend('k_\Omega','k_0^{-2} (E_+ + 0.5 k_\Omega^2 Q_+ - k_\Omega L)^2','Location','southeast')
title('other parameters and E_+ lab frame part of FoM');
xlabel('Iterations');

%%
yend = y(end,:);
kopt = (1/yend(1)) * ( yend(10) + sqrt( yend(10)^2-2*yend(7)*yend(1) ) )

%%
idx = length(z);

fprintf(['Q+ : ',num2str(y(idx,1)),'\n'])
fprintf(['Q- : ',num2str(y(idx,2)),'\n'])
fprintf(['Qx : ',num2str(y(idx,3)),'\n'])
fprintf(['P+ : ',num2str(y(idx,4)),'\n'])
fprintf(['P- : ',num2str(y(idx,5)),'\n'])
fprintf(['Px : ',num2str(y(idx,6)),'\n'])
fprintf(['E+ : ',num2str(y(idx,7)),'\n'])
fprintf(['E- : ',num2str(y(idx,8)),'\n'])
fprintf(['Ex : ',num2str(y(idx,9)),'\n'])
fprintf(['L  : ',num2str(y(idx,10)),'\n'])

%% Plot new FoM4 vs k omega

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);


kk = linspace(-7,-4,1000);
p_opt = a2params(an_h(end,:));
p_noopt = a2params(an_h(1,:));

kk_opt = p_opt(3) / rho;
kk_noopt = p_noopt(3) / rho;

EplusLab = y(end,7) + kk.^2 * y(end,1) / 2 - kk * y(end,10);
EplusLab_opt = y(end,7) + kk_opt.^2 * y(end,1) / 2 - kk_opt * y(end,10);
EplusLab_noopt = y(end,7) + kk_noopt.^2 * y(end,1) / 2 - kk_noopt * y(end,10);

figure; hold on;
plot(kk,EplusLab);
scatter(kk_opt, EplusLab_opt,25,'r');
scatter(kk_noopt, EplusLab_noopt,25,'g');
grid on;

ylabel('FoM term');
xlabel('k_\Omega');
title('(E_+ + k_\Omega^2 Q_+ / 2 - k_\Omega L) vs k_\Omega');

legend('(E_+ + k_\Omega^2 Q_+ / 2 - k_\Omega L) vs k_\Omega','Final optimized k_\Omega','Initial k_\Omega');

global k_perv

kk = linspace(-7,-5,1000);

FoM4 = y(end,7) - kk.^2 * y(end,1) * 0.5 + k_perv;
FoM4_opt = y(end,7) - kk_opt.^2 * y(end,1) / 2 + k_perv;
FoM4_noopt = y(end,7) - kk_noopt.^2 * y(end,1) / 2 + k_perv;

figure; hold on;
plot(kk,FoM4);
scatter(kk_opt, FoM4_opt,25,'r');
scatter(kk_noopt, FoM4_noopt,25,'g');
grid on;

ylabel('FoM term');
xlabel('k_\Omega');
title('(E_+ - k_\Omega^2 Q_+ / 2 + \Lambda) vs k_\Omega , \Lambda = 2.129e-5');

legend('(E_+ - k_\Omega^2 Q_+ / 2 + \Lambda) vs k_\Omega','Final optimized k_\Omega','Initial k_\Omega');


%%

aa = linspace(0.95,1.05,25);
aa2 = linspace(0.75,1.25,25);
FoM_opt = cell(11,1);
FoM_opt_single = cell(11,2);
%aan = an_h(end,:);
%aan = [0.919, 1.554, 0.8655, 1.1153, 0.951, 1.348, 0.893, 1.2827, 1, 1, 1]; 
aan = [1.0, 1.59, 0.95, 1.21, 0.986, 1.34, 0.88, 1.29, 1, 1, 1];
for ii = 1:11
    fprintf([num2str(ii),'/11 \n']);
    atmp = 0;
    if ii == 2
        FoM_opt_tmp = zeros(length(aa2),1);
        for i = 1:length(aa2)
            %fprintf([num2str(i),'/',num2str(length(aa)),'\n']);
            atmp = aan(end,:);
            atmp(ii) = aa2(i) * atmp(ii);
            FoM_opt_tmp(i) = gd_F(atmp);
        end
        FoM_opt{ii} = FoM_opt_tmp;        
    else
        FoM_opt_tmp = zeros(length(aa),1);
        for i = 1:length(aa)
            %fprintf([num2str(i),'/',num2str(length(aa)),'\n']);
            atmp = aan(end,:);
            atmp(ii) = aa(i) * atmp(ii);
            FoM_opt_tmp(i) = gd_F(atmp);
        end
        FoM_opt{ii} = FoM_opt_tmp;        
    end

    aa1 = aan(end,:);
    tmp1 = gd_F(aa1);
    FoM_opt_single{ii,1} = tmp1;
    aa1(ii) = an_h(1,ii);
    tmp1 = gd_F(aa1);
    FoM_opt_single{ii,2} = tmp1;

end
%%
ptitles = {'solenoid position','solenoid strength','quad 1 position','quad 1 strength', ...
    'quad 2 position','quad 2 strength','quad 3 position','quad 3 strength', 'quad 1 rotation', ...
    'quad 2 rotation', 'quad 3 rotation'};

figure;
for ii = 1:11
    subplot(3,4,ii); hold on;
    if ii == 2
        plot(aa2 * an_h(end,ii),FoM_opt{ii});
    else
        plot(aa * an_h(end,ii),FoM_opt{ii});
    end

    scatter(an_h(end,ii), FoM_opt_single{ii,1},25,'r*');
    scatter(an_h(1,ii), FoM_opt_single{ii,2},25,'g');
    grid on;
    title(ptitles{ii});
end

legend('FoM','final','initial');

%%

N = 12;
perChangep = linspace(-0.1,0.1,N);
perChangem = linspace(-0.1,0.1,N);

Fdiff = zeros(N,N);
for i = 1:N
    for j = 1:N
        Fdiff(i,j) = gd_F(an_h(end,:),perChangep(i),perChangem(j));
        fprintf([num2str(i),'|',num2str(j),' / ',num2str(N),'|',num2str(N),'\n']);
    end
end

figure; pcolor(perChangep,perChangem,Fdiff');
colorbar();







