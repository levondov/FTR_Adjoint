function [FoM, FoMp, dQy, dPy, dEy, dLy] = get_F_and_dF(y,i)


global k_solv k_perv k0 k_solvv norm_flag e1 e2 komegaE
global k_solvv
komega = k_solvv;

% calculate adjoint related variables
dPy = zeros(3,length(y(:,1)));
dEy = zeros(3,length(y(:,1)));
dQy = zeros(3,length(y(:,1)));
dLy = zeros(1,length(y(:,1)));

if 0
    % figure of merit broken into pieces for ease of reading version (v paper)
    FoM1 = 0.5*sum(y(i,4:6).^2);
    FoM2 = 0.5*(k0^2)*(y(i,2)^2 + y(i,3)^2);
    FoM3 = 0.5*(k0^(-2))*(y(i,8)^2+y(i,9)^2);
    FoM4 = 0.5*(k0^(-2))*(y(i,7) - 0.5*(komega^2)*y(i,1) + k_perv)^2;
    FoM5 = 0.5*(2*y(i,7)*y(i,1) - y(i,10)^2)^2;    
    
    % adjoint variables calculated from FoM
    dPy(:,i) = y(i,4:6);
    dEy(1,i) = k0^(-2)*(y(i,7)-0.5*komega^2*y(i,1)+k_perv)*0.5*komega^2 - 2*y(i,7)*(2*y(i,7)*y(i,1)-y(i,10)^2);
    dEy(2,i) = -k0^(2)*y(i,2);
    dEy(3,i) = -k0^(2)*y(i,3);
    dQy(1,i) = -k0^(-2)*(y(i,7)-0.5*komega^2*y(i,1)+k_perv) - 2*y(i,1)*(2*y(i,7)*y(i,1)-y(i,10)^2);
    dQy(2,i) = -k0^(-2)*y(i,8);
    dQy(3,i) = -k0^(-2)*y(i,9);
    dLy(:,i) = -2*y(i,10)*(2*y(i,7)*y(i,1)-y(i,10)^2);      
elseif 0
    % figure of merit KISS
    FoM1 = 0.5*(y(i,5)^2 + y(i,6)^2);
    FoM2 = 0.5*(k0^2)*(y(i,2)^2 + y(i,3)^2);
    FoM3 = 0.5*(k0^(-2))*(y(i,8)^2+y(i,9)^2);
    FoM4 = 0.5*(k0^(-2))*e1*( (y(i,7) - 0.5*(komega^2)*y(i,1) + k_perv)^2 + y(i,4)^2 );
    FoM5 = 0.5*e2*( y(i,10) - komega*y(i,1) )^2;

    dPy(1,i) = e1*y(i,4);
    dPy(2,i) = y(i,5);
    dPy(3,i) = y(i,6);
    dEy(1,i) = k0^(-2)*e1*( y(i,7)-0.5*komega^2*y(i,1)+k_perv )*0.5*komega^2 + e2*komega*( y(i,10) - komega*y(i,1) );
    dEy(2,i) = -k0^(2)*y(i,2);
    dEy(3,i) = -k0^(2)*y(i,3);    
    dQy(1,i) = -k0^(-2)*e1*( y(i,7)-0.5*komega^2*y(i,1)+k_perv );
    dQy(2,i) = -k0^(-2)*y(i,8);
    dQy(3,i) = -k0^(-2)*y(i,9);    
    dLy(:,i) = e2*(y(i,10) - komega*y(i,1));
elseif 0
    % figure of merit kiss 2
    % figure of merit KISS
    f5_tmp = ( y(i,7) + 0.5*komega^2*y(i,1) - komega*y(i,10) );
    FoM1 = 0.5*(y(i,5)^2 + y(i,6)^2);
    FoM2 = 0.5*(k0^2)*(y(i,2)^2 + y(i,3)^2);
    FoM3 = 0.5*(k0^(-2))*(y(i,8)^2+y(i,9)^2);
    FoM4 = 0.0; 0.5*(k0^(-2))*e1*( (y(i,7) - 0.5*(komega^2)*y(i,1) + k_perv)^2 + y(i,4)^2 );
    FoM5 = 0.5*e2*(f5_tmp)^2;

    dPy(1,i) = e1*y(i,4);
    dPy(2,i) = y(i,5);
    dPy(3,i) = y(i,6);
    dEy(1,i) = k0^(-2)*e1*( y(i,7)-0.5*komega^2*y(i,1)+k_perv )*0.5*komega^2 + e2*(f5_tmp)*(0.5*komega^2);
    dEy(2,i) = -k0^(2)*y(i,2);
    dEy(3,i) = -k0^(2)*y(i,3);    
    dQy(1,i) = -k0^(-2)*e1*( y(i,7)-0.5*komega^2*y(i,1)+k_perv ) + e2*(f5_tmp)*(-1);
    dQy(2,i) = -k0^(-2)*y(i,8);
    dQy(3,i) = -k0^(-2)*y(i,9);    
    dLy(:,i) = e2*(f5_tmp) * (-1) * (komega);
elseif 1
    % figure of merit kiss 2 no force balance (FoM4)
    f5_tmp = ( y(i,7) + 0.5*komega^2*y(i,1) - komega*y(i,10) );
    f4_tmp = ( y(i,7) - 0.5*(komega^2)*y(i,1) + k_perv );
    FoM1 = 0.5*(y(i,5)^2 + y(i,6)^2) + 0.5*y(i,4)^2;
    FoM2 = 0.5*(k0^2)*(y(i,2)^2 + y(i,3)^2);
    FoM3 = 0.5*(k0^(-2))*(y(i,8)^2+y(i,9)^2);
    FoM4 = 0.5*(k0^(-2))*e1*(f4_tmp)^2;
    FoM5 = (k0^(-2))*0.5*e2*(f5_tmp)^2;

    dPy(1,i) = y(i,4);
    dPy(2,i) = y(i,5);
    dPy(3,i) = y(i,6);
    dEy(1,i) = (k0^(-2))*e2*(f5_tmp)*(0.5*komega^2)*(-1) + (k0^(-2))*e1*(f4_tmp)*(-0.5*komega^2)*(-1);
    dEy(2,i) = -k0^(2)*y(i,2);
    dEy(3,i) = -k0^(2)*y(i,3);    
    dQy(1,i) = (k0^(-2))*e2*(f5_tmp)*(-1) + (k0^(-2))*e1*(f4_tmp)*(-1);
    dQy(2,i) = -k0^(-2)*y(i,8);
    dQy(3,i) = -k0^(-2)*y(i,9);    
    dLy(:,i) = (k0^(-2))*e2*(f5_tmp)*(-1)*(komega);    
end

FoM = (FoM1 + FoM2 + FoM3 + FoM4 + FoM5);
FoMp = [FoM1,FoM2,FoM3,FoM4,FoM5];

end

