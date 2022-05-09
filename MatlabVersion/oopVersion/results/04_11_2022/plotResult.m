figure; hold on;

plot(log10(mo.WY_h),'linewidth',2);
plot(log10(mo3.WY_h),'linewidth',2);
plot(log10(mo1.WY_h),'linewidth',2);
plot(log10(mo2.WY_h),'linewidth',2);


legend('\lambda = 0', ...
    'lambda = 0 , Q_{TE} = 6e-5 , M = 1', ...
    'lambda = 0 , Q_{TE} = 6e-6 , M = 1', ...
    'lambda = 0 , Q_{TE} = 6e-7 , M = 1');

xlabel('Iterations');
grid on;
ylabel('log10(W)');
xlim([0,900]);