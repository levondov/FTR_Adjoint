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






















