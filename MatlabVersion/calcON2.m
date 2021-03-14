function [O,N] = calcON2(z_full, Y, params, k_perv)
%CALCON Summary of this function goes here
global hardedge_flag
% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);

O = cell(length(z_full),1);
N = O;

% iterate through each step and calculate O(z)
for i = 1:length(z_full)
    z = z_full(i);
    
    % solenoid
    z_sol_start = params(1);
    z_sol_end = params(1)+params(2);
    k_sol = 0.0;
    if z >= z_sol_start && z <= z_sol_end % if we are inside solenoid
        k_sol = params(3)/rho; %Ks = Bs / rho = 0.041 (T) / 0.0667
        if ~hardedge_flag
            [xsol,ysol] = solfield(1000,0.025);
            k_sol = interp1(xsol*params(2),ysol*k_sol,z-z_sol_start);
        end  
    end
    
    % quads
    psi = 0;
    k_quad = 0.0;
    z_q1_start = params(4);
    z_q1_end = params(4)+params(5);
    z_q2_start = params(7);
    z_q2_end = params(7)+params(8);
    z_q3_start = params(10);
    z_q3_end = params(10)+params(11);
    
    z_quad_profile = 0:0.05:pi;
    y_quad_profile = cos(z_quad_profile+pi/2).^2;
    scalef = 1; %0.199899;
    scaled = 1;
    scalef_z = params(5); %1mm length
    
    if z >= z_q1_start && z <= z_q1_end % if inside quad 1
        k_quad = params(6)/rho;
        psi = params(13);
        if ~hardedge_flag
            k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q1_start);
        end         
    end
    if z >= z_q2_start && z <= z_q2_end % if inside quad 2
        k_quad = params(9)/rho;
        psi = params(14);
        if ~hardedge_flag
            k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q2_start);
        end         
    end
    if z >= z_q3_start && z <= z_q3_end % if inside quad 3 (== quad 1)
        k_quad = params(12)/rho;
        psi = params(15);
        if ~hardedge_flag
            k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q3_start);
        end         
    end
    
    cq = cos(2*Y(i,11)-2*psi);
    sq = sin(2*Y(i,11)-2*psi);
    
    % space charge stuff
    Q_delta = sqrt( Y(i,1)^2 - Y(i,2)^2 - Y(i,3)^2 );
    ab4 = 1 / Q_delta; % the 4/ab term in equation
    ca_ab4 = -Y(i,2) / ( (Y(i,1)+Q_delta)*Q_delta ); % 4c_alpha/ab
    sa_ab4 = -Y(i,3) / ( (Y(i,1)+Q_delta)*Q_delta ); % 4s_alpha/ab
    
%     % Calculate O and N matrix stuff
%     if (k_sol == 0) && (k_quad == 0)
%         O{i} = [ab4*k_perv, ca_ab4*k_perv, sa_ab4*k_perv;
%             ca_ab4*k_perv, ab4*k_perv, 0;
%             sa_ab4*k_perv, 0, ab4*k_perv];
%         N{i} = [0; -sa_ab4*k_perv; ca_ab4*k_perv];
%     else
%         O{i} = [-k_sol^2/2.0 + ab4*k_perv, 2*k_quad*cq + ca_ab4*k_perv, -2*k_quad*sq + sa_ab4*k_perv;
%             2*k_quad*cq + ca_ab4*k_perv, -k_sol^2/2.0 + ab4*k_perv, 0;
%             -2*k_quad*sq + sa_ab4*k_perv, 0, -k_sol^2/2.0 + ab4*k_perv];
% 
%         N{i} = [0; 2*k_quad*sq - sa_ab4*k_perv; 2*k_quad*cq + ca_ab4*k_perv];
%     end
        O{i} = [-k_sol^2/2.0 + ab4*k_perv, 2*k_quad*cq + ca_ab4*k_perv, -2*k_quad*sq + sa_ab4*k_perv;
            2*k_quad*cq + ca_ab4*k_perv, -k_sol^2/2.0 + ab4*k_perv, 0;
            -2*k_quad*sq + sa_ab4*k_perv, 0, -k_sol^2/2.0 + ab4*k_perv];

        N{i} = [0; 2*k_quad*sq - sa_ab4*k_perv; 2*k_quad*cq + ca_ab4*k_perv];
    
    
    
end

end

