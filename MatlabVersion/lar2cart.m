function [y] = lar2cart(dYdt2,dphi)
%LAR2CART Summary of this function goes here

% Larmor to fixed frame transform
lc = cos(2*dYdt2(:,11));
ls = sin(2*dYdt2(:,11));

N = length(dYdt2(:,1));
dYdt = zeros(N,11);

dYdt(:,1) = dYdt2(:,1);
dYdt(:,2) = dYdt2(:,2).*lc - dYdt2(:,3).*ls;
dYdt(:,3) = dYdt2(:,2).*ls + dYdt2(:,3).*lc;
dYdt(:,4) = dYdt2(:,4);
dYdt(:,5) = dYdt2(:,5).*lc - dYdt2(:,6).*ls + 2*dphi.*( -dYdt2(:,2).*ls - dYdt2(:,3).*lc );
dYdt(:,6) = dYdt2(:,5).*ls + dYdt2(:,6).*lc + 2*dphi.*( dYdt2(:,2).*lc - dYdt2(:,3).*ls );
dYdt(:,7) = dYdt2(:,7) + 2*(dphi.^2).*dYdt2(:,1) + 2*dphi.*dYdt2(:,10);
dYdt(:,8) = ( dYdt2(:,8).*lc - dYdt2(:,9).*ls ) + 2*dphi.*( -dYdt2(:,5).*ls - dYdt2(:,6).*lc  ) + 2*(dphi.^2).*( -dYdt2(:,2).*lc + dYdt2(:,3).*ls );
dYdt(:,9) = ( dYdt2(:,8).*ls + dYdt2(:,9).*lc ) + 2*dphi.*( dYdt2(:,5).*lc - dYdt2(:,6).*ls  ) + 2*(dphi.^2).*( -dYdt2(:,2).*ls - dYdt2(:,3).*lc );
dYdt(:,10) = dYdt2(:,10) + 2*dphi.*dYdt2(:,1);
dYdt(:,11) = dYdt2(:,11);

y=dYdt;



end

