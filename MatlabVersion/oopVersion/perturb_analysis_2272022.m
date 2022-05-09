%% Initial conditions

% 20mA
X0 = [
   0.007363052819400*(1-0.27)
  -0.003144679253067*(1-0.25)
                   0
   0.000000302378449
   0.000000350230487
                   0
   0.554383041347526
   0.239127669593029
                   0
                   0
                   0] * 1e-4;      

% 10 mA
% X0 = [
%    0.035153540129256
%   -0.004652095581008
%                    0
%    0.000000000011345
%    0.000000000011576
%                    0
%    0.167236580717593
%    0.022131499566718
%                    0
%                    0
%                    0] * 1e-4;   
                              
% 5 mA
% X0 = [
%    0.023821680783865
%   -0.003156661821524
%                    0
%    0.000000000011829
%    0.000000000014758
%                    0
%    0.246801907248502
%    0.032704248082207
%                    0
%                    0
%                    0] * 1e-4;                  

% 0mA
% X0 = [
%    0.015000620213669
%   -0.001995219530101
%                    0
%    0.000000000011346
%    0.000000000016410
%                    0
%    0.391985530397016
%    0.052137656618415
%                    0
%                    0
%                    0] * 1e-4;   
               
% calculate E+ and E- for a give emittance
ex2 = 6.2e-6^2; % target emittance X
ey2 = 6.2e-6^2; % target emittance Y
X0 = CalcValsForEmittance(X0,ex2,ey2);

% create lattice
turns = 1024;
sidx = 1;
eidx = 1024;
plotidxstart = 2;
plotidxend = 1024;

% setup moment object
mom = MomentSolverPeriodic(10e3, 20e-3, X0); % energy, beam current, initial conditions
mom.h = 100; % integration steps in each element
% create lattice
an = ones(5,1)';
mom = CreateLattice(mom, an, turns, 0);

% run moment equations
mom = mom.RunMoments(true);

%% Extract out beam size at center of quadrupole

% grab quad center locations
numQuads = turns*2; % 2 quad per a period
quad2quadspacing = 0.04;
[quadlocIdx, quadlocs] = FindQuadCenters(mom, quad2quadspacing, numQuads);

%% Freq analysis
cspeed = 299792458;
% calculate beam size (<x^2> <y^2>) at quad centers
% every other quad only
zb = quadlocs(1:2:end);
tb = zb/cspeed;
xb = mom.y(quadlocIdx(1:2:end),1) + mom.y(quadlocIdx(1:2:end),2);
yb = mom.y(quadlocIdx(1:2:end),2) - mom.y(quadlocIdx(1:2:end),2);

% verify in plot
if 0
    % red line across black lattice
    mom.PlotBeamSize();
    for i = 1:length(quadlocs)
        plot([quadlocs(i),quadlocs(i)],[-0.5e-6,0.5e-6],'r--')
    end
    % mark every peak/trough in beamsize with red dot
    plot(zb,xb,'r.');
    plot(zb,yb,'r.');    
    legend('off');
end

% take an fft of the data
[xb_freq,xb_mag] = freqAnalysis(zb(sidx:eidx),xb(sidx:eidx));
[yb_freq,yb_mag] = freqAnalysis(zb(sidx:eidx),yb(sidx:eidx));

%% Perturbed case

%Perturb initial conditions
X0p = X0;
X0p(3) = X0(3) * 1.10; % 5 percent

% update moment solver initial conditions and rerun
% also make a copy
mom_p = mom;
mom_p.initialMoments = X0p;
mom_p = mom_p.RunMoments(true);

%% Freq analysis pert

% calculate beam size (<x^2> <y^2>) at quad centers
% look at every other quad
zb_p = quadlocs(1:2:end);
tb_p = zb/cspeed;
xb_p = mom_p.y(quadlocIdx(1:2:end),1) + mom_p.y(quadlocIdx(1:2:end),2);
yb_p = mom_p.y(quadlocIdx(1:2:end),1) - mom_p.y(quadlocIdx(1:2:end),2);

% take an fft of the data
[xb_p_freq,xb_p_mag] = freqAnalysis(zb_p(sidx:eidx),xb_p(sidx:eidx));
[yb_p_freq,yb_p_mag] = freqAnalysis(zb_p(sidx:eidx),yb_p(sidx:eidx));

%% Plot the results
figure; hold on;
%subplot(3,1,3); hold on;

NN = length(zb);

plot(xb_freq(plotidxstart:plotidxend),xb_mag(plotidxstart:plotidxend),'.-');
plot(yb_freq(plotidxstart:plotidxend),yb_mag(plotidxstart:plotidxend),'^-');
plot(xb_p_freq(plotidxstart:plotidxend),xb_p_mag(plotidxstart:plotidxend),'>-');
plot(yb_p_freq(plotidxstart:plotidxend),yb_p_mag(plotidxstart:plotidxend),'s-');

% plot results
legend('x^2','y^2','x^2 pert','y^2 pert');
xlabel('k wavenumber [rad/meters]'); ylabel('Mag [a.u.]');
title('20 mA, Q+ Pert');



%% Functions
function [zlocIdx,zloc] = FindQuadCenters(mom,q2qSpacing, totalQuads)

zlocIdx = [];
zloc = [];
quadCentersZ = 0:q2qSpacing:(totalQuads*q2qSpacing);
N = length(quadCentersZ);

% loop through each quad locations
for i = 1:N
        % find closest index value to that center
        idx = find(mom.z >= quadCentersZ(i),1);
        if ( ~isempty(idx) )
            zlocIdx(end+1) = idx;
            zloc(end+1) = mom.z(idx);
        end
end
end

function [X] = CalcValsForEmittance(X0,ex2,ey2)

% calculate E+ and E- for a give emittance
% ex2 = 7.6e-6^2; % target emittance X
% ey2 = 7.6e-6^2; % target emittance Y

denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
X0(7) = Ep;
X0(8) = Em;
X = X0';
end



               