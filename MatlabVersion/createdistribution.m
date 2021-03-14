function [x0,y0,xp0,yp0] = createdistribution(Np)
% gaussian distribution in x and y
xmean = 0.0;
ymean = 0.0;
xsigma = 2.26e-3;
ysigma = 0.226e-3;
x0 = normrnd(xmean,xsigma,[Np,1]);
y0 = normrnd(ymean,ysigma,[Np,1]);
% random distr in xp,yp
xpmin = -0.0;
xpmax = 0.0;
ypmin = -0.0;
ypmax = 0.0;
%xp0 = ((xpmax-xpmin)*rand([Np,1])+xpmin)*0.0;
%yp0 = ((ypmax-ypmin)*rand([Np,1])+ypmin)*0.0;
xp0 = normrnd(xmean,0.013,[Np,1]);
yp0 = normrnd(ymean,0.013,[Np,1]);

% manually manipulate  data
if 0
    th = zeros(Np,1); rr=zeros(Np,1);
    tt_pert = 0.01;
    for ii = 1:Np
        [th(ii),rr(ii)] = cart2pol(x0(ii),y0(ii));
        th(ii) = th(ii) + tt_pert;
        [x2,y2] = pol2cart(th(ii),rr(ii));
        xl = x2-x0(ii);
        yl = y2-y0(ii);     
        xp0(ii) = 1e-3*xl*1e6;
        yp0(ii) = 1e-3*yl*1e6;  
    end
end

end