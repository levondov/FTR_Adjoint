   function yout = ode2(F,t0,h,tfinal,y0)
   % ODE4  Classical Runge-Kutta ODE solver.
   %   yout = ODE4(F,t0,h,tfinal,y0) uses the classical
   %   Runge-Kutta method with fixed step size h on the interval
   %      t0 <= t <= tfinal
   %   to solve
   %      dy/dt = F(t,y)
   %   with y(t0) = y0.

   %   Copyright 2014 - 2015 The MathWorks, Inc.
   
      y = y0;
      yout = y;
      for t = t0 : h : tfinal-h
         s1 = F(t,y);
         s2 = F(t+2*h/3, y+2*h*s1./3);
         y = y + h*(s1/4 + 3*s2/4);
         yout = [yout, y]; %#ok<AGROW>
      end
      
      yout=yout';
