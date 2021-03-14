function [dYdt] = odefcn(z,Y,params)
%
% Y = [Q,P,E,L]
% Y = [Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L, phi]
%     [   1       2    3     4       5     6    7      8     9  10 11

% Parameters
e         = 1.60217733E-19; %C
m         = 9.1093897E-31; %kg
Energy    = 5e3;%9967.08077; % eV
c         = 2.997924E8; % m/s

gamma     = 1+((Energy)/(510998.9461));
beta      = sqrt((gamma*gamma)-1)/gamma;
bg           = beta*gamma;
rho         = bg*c*(m/e);

% hardedge or not
global hardedge_flag

% just diagnostics to grab more data
global extra_params_flag extra_params_flag_firstcall

% flag for space charge
global k_perv k_solvv

if extra_params_flag
    global zv cqv sqv k_quadv k_solv O N
    if extra_params_flag_firstcall
        % first time calling the function
        zv = []; cqv = []; sqv = []; k_quadv = []; k_solv = [];O = {}; N = {};
        extra_params_flag_firstcall = 0;
    end
end

%% important constants for calculation
% params = [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
% strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]

% solenoid
z_sol_start = params(1);
z_sol_end = params(1)+params(2);
k_sol = 0.0;
k_solvv = params(3)/rho; % just need max solenoid value for GD optimizations
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
%% Latice ends at z = 5.564
% non hardedge field stuff
z_quad_profile = 0:0.05:pi;
y_quad_profile = cos(z_quad_profile+pi/2).^2;
scalef = 1; %0.296097; %0.199899;
scaled = 1; %0.296097;
scalef_z = params(5); %1mm length
%%%%%%%%%
if z >= z_q1_start && z <= z_q1_end % if inside quad 1
    k_quad = params(6)/rho;
    psi = params(13);
    if ~hardedge_flag
        k_quad = k_quad*scalef;
        k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q1_start);
    end
end
if z >= z_q2_start && z <= z_q2_end % if inside quad 2
    k_quad = params(9)/rho;
    psi = params(14);
    if ~hardedge_flag
        k_quad = k_quad*scalef;
        k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q2_start);
    end
end
if z >= z_q3_start && z <= z_q3_end % if inside quad 3 (== quad 1)
    k_quad = params(12)/rho;
    psi = params(15);
    if ~hardedge_flag
        k_quad = k_quad*scalef;
        k_quad = interp1([-1e-6,linspace(0,1,length(z_quad_profile)),1.0+1e-6]*scalef_z,[0,y_quad_profile*k_quad,0],z-z_q3_start);
    end
end
%psi=0;
cq = cos(2*Y(11)-2*psi);
sq = sin(2*Y(11)-2*psi);

% space charge stuff
Q_delta = sqrt( Y(1)^2 - Y(2)^2 - Y(3)^2 );
ab4 = 1 / Q_delta; % the 4/ab term in equation
ca_ab4 = -Y(2) / ( (Y(1)+Q_delta)*Q_delta ); % 4c_alpha/ab
sa_ab4 = -Y(3) / ( (Y(1)+Q_delta)*Q_delta ); % 4s_alpha/ab

% Calculate O and N matrix stuff
O_mat = [-k_sol^2/2.0 + ab4*k_perv, 2*k_quad*cq + ca_ab4*k_perv, -2*k_quad*sq + sa_ab4*k_perv;
    2*k_quad*cq + ca_ab4*k_perv, -k_sol^2/2.0 + ab4*k_perv, 0;
    -2*k_quad*sq + sa_ab4*k_perv, 0, -k_sol^2/2.0 + ab4*k_perv];

N_mat = [0; 2*k_quad*sq - sa_ab4*k_perv; 2*k_quad*cq + ca_ab4*k_perv];

if extra_params_flag
    zv(end+1) = z;
    cqv(end+1) = cq;
    sqv(end+1) = sq;
    k_quadv(end+1) = k_quad;
    k_solv(end+1) = k_sol;
    O{end+1} = O_mat;
    N{end+1} = N_mat;
end


%% System of 10 equations to solve

% dQ/dz
dY1dt = [ Y(4);
    Y(5);
    Y(6) ];

% dP/dz
dY2dt = Y(7:9) + O_mat*Y(1:3);

% dE/dz
dY3dt = O_mat*Y(4:6) + N_mat*Y(10);

% dL/dz
dY4dt = -N_mat'*Y(1:3);

% d phi / dz
dY5dt = -1*k_sol/2.0;

% put all together
dYdt2 = [dY1dt; dY2dt; dY3dt; dY4dt; dY5dt];
dYdt = dYdt2;

end









