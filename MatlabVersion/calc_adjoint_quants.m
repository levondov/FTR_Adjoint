function [Mq,Mp,Mn] = calc_adjoint_quants(Y,k_perv)
%CALC_ADJOINT_QUANTS Summary of this function goes here

%Y(1) - Q+
%Y(2) - Q-
%Y(3) - Qx
%Y(4) - P+
%Y(5) - P-
%Y(6) - Px
%Y(7) - E+
%Y(8) - E-
%Y(9) - Ex
%Y(10) - L


% calculate all the required quantities for adjoint equations
Q_delta = sqrt( Y(1)^2 - Y(2)^2 - Y(3)^2 );
H = k_perv*log(Y(1) + Q_delta);
Q_deltaplus = Q_delta + Y(1);

V1 = -(k_perv/(Q_delta^2))*[Y(4)-Y(2)*Y(5)/Q_deltaplus-Y(3)*Y(6)/Q_deltaplus,
    -Y(2)*Y(4)/Q_deltaplus+Y(5),
    -Y(3)*Y(4)/Q_deltaplus+Y(6)];
U1t = (1/Q_delta)*[Y(1), -Y(2), -Y(3)];

V2 = (k_perv/(Q_delta*(Q_deltaplus^2)))*[Y(2)*Y(5)+Y(3)*Y(6),Y(2)*Y(4),Y(3)*Y(4)]';
U2t = U1t + [1,0,0];

V3 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(5), Y(4), 0]';
V4 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(6), 0, Y(4)]';
U3t = [0, 1, 0];
U4t = [0,0,1];

Mp = V1*U1t + V2*U2t + V3*U3t + V4*U4t;
%%%%

W1 = -(k_perv/Q_delta)*[1, Y(2)/Q_deltaplus, Y(3)/Q_deltaplus]';
W2 = (k_perv/(Q_delta*(Q_deltaplus^2)))*[Y(2)^2+Y(3)^2,Y(2)*Y(1),Y(3)*Y(1)]';
W3 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(2),Y(1),0]';
W4 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(3),0,Y(1)]';

Mq = W1*U1t + W2*U2t + W3*U3t + W4*U4t;
%%%%

X1 = -k_perv*[0, Y(3)/(Q_deltaplus*(Q_delta^2)), -Y(2)/(Q_deltaplus*(Q_delta^2))]';
X2 = -k_perv*[0, Y(3)/(Q_delta*(Q_deltaplus^2)), -Y(2)/(Q_delta*(Q_deltaplus^2))]';
X3 = k_perv*[0,0,-1/(Q_delta*Q_deltaplus)]';
X4 = k_perv*[0,1/(Q_delta*Q_deltaplus),0]';

Mn = X1*U1t + X2*U2t + X3*U3t + X4*U4t;


end









