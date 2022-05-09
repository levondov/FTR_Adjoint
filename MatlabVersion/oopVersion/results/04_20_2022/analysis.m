%%
%
%
%
%


%% First run Opt like normal.
%

mom.PlotBeamSize;
mom.PlotBeamEmit;

X0 = [
       0.000407212721295
  -0.000183484869158
                   0
  -0.000000002299260
  -0.000000001843467
                   0
   0.118445744827765
   0.053370145026373
                   0
                   0
                   0] * 1e-3;

%% 
% Next we Find QTE

L = 0.08 * 2; % 2 periods
M = 2;
Qp = mom.y(:,1);
z = mom.z;

QTE = linspace(6e-8,6e-6,500);
F = zeros(length(QTE),1);

for i = 1:length(QTE)
    fs = (1/L) * (Qp./QTE(i) - 1).^M;
    F(i) = trapz(z,fs);
end

figure;
plot(QTE,log10(F),'linewidth',2);
xlabel('QTE'); ylabel('log10(FE)'); title('log10(FE) vs QTE'); grid on;

[v,idx] = min(F);
QTEbest = QTE(idx);


