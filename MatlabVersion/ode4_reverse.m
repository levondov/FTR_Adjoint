function [yout2,yout1] = ode4_reverse(F1,F2,h,y01,y02,t0,tfinal)
% ODE4  Classical Runge-Kutta ODE solver.
%   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
%   Runge-Kutta method with fixed step size h on the interval
%      t0 <= t <= tfinal
%   to solve
%      dy/dt = F(t,y)
%   with y(t0) = y0.

%   Copyright 2014 - 2015 The MathWorks, Inc.

y1 = y01;
y2 = y02;
yout1 = y1;
yout2 = y2;

% first do non-adjoint integration backwards at h/2
for t = t0 : h/2.0 : tfinal-h/2.0
    s1 = F1(t, y1);
    s2 = F1(t+h/2, y1+h.*s1./2);
    s3 = F1(t+h/2, y1+h.*s2./2);
    s4 = F1(t+h, y1+h.*s3);
    y1 = y1 + h.*(s1 + 2*s2 + 2*s3 + s4)./6;
    yout1 = [yout1, y1]; %#ok<AGROW>
end

% next do adjoint integration backwards
idx_loc = 1;
for t = t0 : h : tfinal-h
    s1 = F2(t,y2, yout1(:,idx_loc));
    s2 = F2(t+h/2, y2+h.*s1./2, yout1(:,idx_loc+1));
    s3 = F2(t+h/2, y2+h.*s2./2, yout1(:,idx_loc+1));
    s4 = F2(t+h, y2+h.*s3, yout1(:,idx_loc+1));
    y2 = y2 + h.*(s1 + 2*s2 + 2*s3 + s4)./6;
    yout2 = [yout2, y2]; %#ok<AGROW>
    idx_loc = idx_loc + 2;
end

yout2=yout2';
yout1=yout1';


end
