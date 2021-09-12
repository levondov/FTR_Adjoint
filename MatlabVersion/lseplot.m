function [ output_args ] = lseplot(z,y,motion,tstring)

%% Now plot everything to compare
figure('units','pixels','position',[100,100,750,800]); 
subplot(3,2,1); hold on;
plot(z,y(:,1),'k','linewidth',3);
plot(z,y(:,2),'b','linewidth',3);
plot(z,y(:,3),'r','linewidth',3);
legend('q+','q-','qx','location','northwest');
%ylim([-5e-5,5e-5]);
ylabel('Q+,Q-,Qx [m^2]'); xlabel('z [m]'); grid on; ylim([-7e-6,7e-6]);

subplot(3,2,2); hold on;
plot(z,y(:,4),'k','linewidth',3);
plot(z,y(:,5),'b','linewidth',3);
plot(z,y(:,6),'r','linewidth',3);
legend('p+','p-','px','location','northwest');
%ylim([-5e-5,5e-5]);
ylabel('P+,P-,Px [m]'); xlabel('z [m]'); grid on; ylim([-1e-4,1e-4]);

global k_solv zv
kidx = find(k_solv ~= 0,1);
zidx = find(z >= zv(kidx),1) - 1;
komega = zeros(length(z),1);
komega(zidx:end) = k_solv(kidx);

subplot(3,2,3); hold on;
plot(z,y(:,7),'k','linewidth',3);
plot(z,y(:,8),'b','linewidth',3);
plot(z,y(:,9),'r','linewidth',3);
plot(z,y(:,7) + 0.5*komega.^2.*y(:,1) - komega.*y(:,10),'g--','linewidth',3);

legend('e+','e-','ex','e+_{lab}','location','northwest');
%ylim([-2e-5,3.5e-5]);
ylabel('E+,E-,Ex '); xlabel('z [m]'); grid on; ylim([-6e-4,6e-4]);

subplot(3,2,4); hold on;
plot(z,y(:,10),'k','linewidth',3);
plot(z,motion,'b','linewidth',3);
plot(z,y(:,10)-komega.*y(:,1),'g','linewidth',3,'linestyle','-.');

%plot(z,y(:,11),'r','linewidth',3);
legend('L','C','L_{lab}','location','northwest');
%ylim([-3e-5,5e-6]);
ylabel('L,constant of motion'); xlabel('z [m]'); grid on; ylim([-6e-5,2e-5]);

% figure('units','pixels','position',[100,100,750,500]); 
% subplot(2,2,1); hold on;
% plot(z,dPy(1,:),'k','linewidth',3);
% plot(z,dPy(2,:),'b','linewidth',3);
% plot(z,dPy(3,:),'r','linewidth',3);
% legend('\delta Py +','\delta Py - ','\delta Py x','location','northwest');
% %ylim([-5e-5,5e-5]);
% ylabel('\delta Py'); xlabel('z [m]');
% 
% subplot(2,2,2); hold on;
% plot(z,dEy(1,:),'k','linewidth',3);
% plot(z,dEy(2,:),'b','linewidth',3);
% plot(z,dEy(3,:),'r','linewidth',3);
% legend('\delta Ey +','\delta Ey -','\delta Ey x','location','northwest');
% %ylim([-5e-5,5e-5]);
% ylabel('\delta Ey'); xlabel('z [m]');
% 
% subplot(2,2,3); hold on;
% plot(z,dQy(1,:),'k','linewidth',3);
% plot(z,dQy(1,:),'b','linewidth',3);
% plot(z,dQy(1,:),'r','linewidth',3);
% legend('\delta Qy +','\delta Qy -','\delta Qy x','location','northwest');
% %ylim([-5e-5,5e-5]);
% ylabel('\delta Qy'); xlabel('z [m]');
% 
% subplot(2,2,4); hold on;
% plot(z,dLy(1,:),'k','linewidth',3);
% legend('\delta Ly','location','northwest');
% %ylim([-5e-5,5e-5]);
% ylabel('\delta Ly'); xlabel('z [m]');
global zv cqv sqv k_quadv k_solv k_perv O N

subplot(3,2,5); hold on;
plot(z,y(:,11),'r','linewidth',3);
plot(zv,cqv,'-','Linewidth',1); hold on; 
plot(zv,sqv,'-','Linewidth',1);
legend('\phi','c_q','s_q','location','best');
%ylim([-3e-5,5e-6]);
ylabel('\phi , c_q, s_q'); xlabel('z [m]');
title('c_q, s_q, and \phi')


subplot(3,2,6); hold on;
plot(zv,k_quadv/1e2,'k-','Linewidth',1);
plot(zv,k_solv,'m-','Linewidth',1);
%xlim([0,z_interval(2)])
title('Magnet profiles');
legend('k_{quad}/10','k_{sol}','location','best');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', tstring, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize',20)

figure; hold on;
plot(z,y(:,1)+y(:,2));
plot(z,y(:,1)-y(:,2));

plot(zv,k_quadv/1e2*1e-6*0.2,'k-','Linewidth',1);
plot(zv,k_solv*1e-6*0.25,'m-','Linewidth',1);

xlabel('Z position (m)'); ylabel('Moments [m^2]'); title('Beam size through FTR transformer');
legend('\langle x^2 \rangle','\langle y^2 \rangle','K_{quad}','k_{solenoid}');
grid on;

end

