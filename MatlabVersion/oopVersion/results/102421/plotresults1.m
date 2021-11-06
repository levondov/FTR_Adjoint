
files = {
    'matching_0mA_aperiodic_5perc.mat'
    'matching_1mA_aperiodic_5perc.mat'
    'matching_2mA_aperiodic_5perc.mat'
    'matching_5mA_aperiodic_5perc.mat'
    'matching_10mA_aperiodic_5perc.mat'
    'matching_20mA_aperiodic_5perc.mat'
    'matching_40mA_aperiodic_5perc.mat'
    'matching_60mA_aperiodic_5perc.mat'
    'matching_80mA_aperiodic_5perc.mat'
    'matching_825mA_aperiodic_5perc.mat'
    'matching_85mA_aperiodic_5perc.mat'
    'matching_100mA_aperiodic_5perc.mat'
};

N = length(files);

xs = zeros(N,2);
ys = zeros(N,2);

for i = 1:N
    fprintf([num2str(i),'/',num2str(N),'\n']);
   load(files{i}); 
   mom = CreateLatticeAperiodic(mom, an, 5, 5, 1);
   xi = mom.y(1,1) + mom.y(1,2);
   xf = mom.y(end,1) + mom.y(end,2);
   yi = mom.y(1,1) - mom.y(1,2);
   yf = mom.y(end,1) - mom.y(end,2);   
   xs(i,:) = [xi,xf];
   ys(i,:) = [yi,yf];
end

%%

crnts = [0,1,2,5,10,20,40,60,80,82.5,85,100];

figure; hold on;

plot(crnts,log10((abs(xs(:,2)-xs(:,1)))./xs(:,1)),'.-');
plot(crnts,log10((abs(ys(:,2)-ys(:,1)))./ys(:,1)),'.-');

ylabel('log10[(xf^2-xi^2)/xi^2]');
xlabel('Beam current [mA]');
grid on;