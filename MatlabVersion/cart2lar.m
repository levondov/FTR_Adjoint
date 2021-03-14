function [y_lar] = cart2lar(y,dphi)
%CART2LAR Summary of this function goes here
%   convert from cartesian to larmor
% y should be Nx11 in cartesian form

N = length(y(:,1));
y_lar = zeros(N,11);

lc = cos(2*y(:,11));
ls = sin(2*y(:,11));

% Q+,Q-,Qx
y_lar(:,1) = y(:,1);
y_lar(:,2) = y(:,2).*lc + y(:,3).*ls;
y_lar(:,3) = y(:,3).*lc - y(:,2).*ls;

% P+,P-,Px
y_lar(:,4) = y(:,4);
y_lar(:,5) = y(:,5).*lc + y(:,6).*ls + 2*dphi.*y_lar(:,3);
y_lar(:,6) = -y(:,5).*ls + y(:,6).*lc - 2*dphi.*y_lar(:,2);

% L
y_lar(:,10) = y(:,10) - 2*dphi.*y_lar(:,1);

% E+,E-,Ex
d = -y_lar(:,5).*ls - y_lar(:,6).*lc;
e = -y_lar(:,2).*lc + y_lar(:,3).*ls;
f = y_lar(:,5).*lc - y_lar(:,6).*ls;
g = -y_lar(:,2).*ls - y_lar(:,3).*lc;

y_lar(:,7) = y(:,7) - 2*(dphi.^2).*y_lar(:,1) - 2*dphi.*y_lar(:,10);
y_lar(:,8) = y(:,8).*lc + y(:,9).*ls + 2.*dphi.*y_lar(:,6) + 2.*(dphi.^2).*y_lar(:,2);
y_lar(:,9) = -y(:,8).*ls + y(:,9).*lc - 2.*dphi.*y_lar(:,5) + 2.*(dphi.^2).*y_lar(:,3);

% phi
y_lar(:,11) = y(:,11);


end

