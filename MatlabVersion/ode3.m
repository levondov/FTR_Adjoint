   function yout = ode3(F,t0,h,tfinal,y0)
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
         s1 = h.*F(t,y);
         s2 = h.*F(t+h/2, y+s1./2);
         s3 = h.*F(t+h, y-s1+2*s2);
         y = y + (s1 + 4*s2 + s3)./6;
         yout = [yout, y]; %#ok<AGROW>
      end
      
      yout=yout';
