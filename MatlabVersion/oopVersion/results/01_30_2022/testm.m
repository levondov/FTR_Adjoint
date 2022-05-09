[NN,N] = size(mom.y);

fomh = zeros(NN,1);
yy = mom.y;

for ii = 1:NN
    mom.initialMoments = yy(ii,:)';
    mom = mom.RunMoments();
   [f0,f0p] = mom.GetFAndDF1();
   fomh(ii) = f0;
   fprintf([num2str(ii),'/',num2str(NN),'\n']);
end

%%


files = {
    '00_11mA_jcorr_em76.mat',
    '00_12mA_jcorr_em76.mat',
    '00_13mA_jcorr_em76.mat',
    '00_135mA_jcorr_em76.mat',
    '00_14mA_jcorr_em76.mat'    
};
ln = {'11mA','12mA','13mA','13.5mA','14mA'};
nn = {'Q+','Q-','Qx','P+','P-','Px','E+','E-','Ex'};
NN = 5;

figure; hold on;
for i = 1:NN
    load(files{i}); 
    
    for j = 1:9
        subplot(3,3,j); hold on;        
        plot(Xn_h(:,j),'linewidth',2);
        title(nn{j});
    end
end

subplot(3,3,1);
legend(ln);

figure; hold on;
for i = 1:NN
    load(files{i}); 
    
    plot(log10(f_h),'linewidth',2);
end
legend(ln);
title('FoM vs iterations');

figure; hold on;
for i = 1:NN
    load(files{i}); 
    
    plot(j_h(:,1),'linewidth',2);
end
legend(ln);
title('J vs iterations');

figure; hold on;
for i = 1:NN
    load(files{i}); 
    
    for j = 1:9
        subplot(3,3,j); hold on;        
        plot(df_h(j,:),'linewidth',2);
        title(nn{j});
    end
end

subplot(3,3,1);
legend(ln);
title('dF/dQ+ vs iterations');

figure; hold on;
for i = 1:NN
    load(files{i}); 
    
    for j = 1:9
        tmp1 = zeros(length(dj_h),1);
        for jj = 1:length(tmp1)
            tmp1(jj) = dj_h{jj}(j,1);
        end        
        subplot(3,3,j); hold on;        
        plot(tmp1,'linewidth',2);
        title(nn{j});
    end
end

subplot(3,3,1);
legend(ln);
title('dJ/dQ+ vs iterations');






