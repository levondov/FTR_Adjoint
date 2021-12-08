classdef MomentSolverPeriodic
    %MOMENTSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physics parameters
        energy % [eV]
        current % [A]
        pipeRadius % [m]
        rigidity
        perveance
        
        % lattice variables
        lattice
        zstart
        zend
        
        % integration variables
        h = 10000; % step size
        
        % Moment variables
        initialMoments
        z
        y
        k
        com
        yOri
        yAdj
        zAdj
        kAdj
        
        % Moment equations
        k0 = 7;
        Omat
        Nmat
        
    end
    
    methods
        function obj = MomentSolverPeriodic(energy, current, initialMoments)
            obj.energy = energy;
            obj.current = current;
            obj.initialMoments = initialMoments;
            [obj.perveance, obj.rigidity] = obj.CalculateBeamProperties(obj.energy, obj.current);
        end
        
        
        function [perv, rho] = CalculateBeamProperties(obj, energy, current)
            % Parameters
            e         = 1.60217733E-19; %C
            m         = 9.1093897E-31; %kg
            Energy    = energy;%9967.08077; % eV
            c         = 2.997924E8; % m/s
            
            gamma     = 1+((Energy)/(510998.9461));
            beta      = sqrt((gamma*gamma)-1)/gamma;
            bg        = beta*gamma;
            rho       = bg*c*(m/e);
            v         = beta*c;
            
            perv = (1/(4*pi))*(c*377*current) / (m*v^3*gamma^3/e);
        end
        
        % Create the accelerator lattice
        function obj = CreateLatticeProfile(obj, dB, qstart, qend, qrot, zstart, zend, repeat, verbose)
            if nargin == 7
                verbose = false;
            end
            
            % setup lattice
            lattice = [];
            % drift space before first quad
            if ( (qstart(1) - zstart) > 0 )
                lattice(end+1,:) = [zstart, qstart(1), 0.0, 0.0];
            end
            
            for i = 1:length(dB)
                % add quad
                lattice(end+1,:) = [qstart(i), qend(i), dB(i), qrot(i)];
                % add drift to next quad
                if (i < length(dB))
                    lattice(end+1,:) = [qend(i), qstart(i+1), 0.0, 0.0];
                end
            end
            % add final drift
            if ( (zend - qend(end)) > 0 )
                lattice(end+1,:) = [qend(end), zend, 0.0, 0.0];
            end
            
            % repeating lattice elements
            latticeCopy = lattice;
            for i = 1:(repeat-1)
                latticeCopy(:,[1,2]) = latticeCopy(:,[1,2]) + zend;
                lattice = [lattice; latticeCopy];
            end
            
            % assign properties
            obj.lattice = lattice;
            obj.zstart = zstart;
            obj.zend = zend * repeat;
            
            if verbose
                lattice
            end
        end
        
        function obj = RunMoments(obj, verbose)
            if nargin < 2
                verbose = false;
            end
            
            %[z,y,kval,Omat,Nmat] = obj.ode3(@(t,Y) obj.odefcn(t,Y), obj.zstart, obj.h, obj.zend, obj.initialMoments, verbose);
            [z,y,kval,Omat,Nmat] = obj.ode3(@(t,Y,db,rot) obj.OdeMoments(t,Y,db,rot), obj.h, obj.initialMoments, verbose);
            
            obj.z = z;
            obj.y = y;
            obj.k = kval;
            obj.Omat = Omat;
            obj.Nmat = Nmat;
            
            % Constant of Motion 1% 0.5*Tr(J_4^2 sigma^2)
            L = [y(:,10)];
            EQ = y(:,7).*y(:,1) + y(:,8).*y(:,2) + y(:,9).*y(:,3);
            PP = y(:,4).^2 + y(:,5).^2 + y(:,6).^2;
            motion1 = EQ + (1/2)*L.^2 - (1/2)*PP;            
            obj.com = motion1;
        end
        
        function obj = RunMomentsAdjoint(obj, verbose)
            if nargin < 2
                verbose = false;
            end
            
            % grab adjoint gradients
            [~,~,dQy, dPy, dEy, dLy] = obj.GetFAndDF1();
            initialMoments = [dQy; dPy; dEy; dLy; 0.0; obj.y(end,:)'];
            
            %[z,y,kval,Omat,Nmat] = obj.ode3(@(t,Y) obj.odefcn(t,Y), obj.zstart, obj.h, obj.zend, obj.initialMoments, verbose);
            [z,y,kval,Omat,Nmat] = obj.ode3(@(t,Y,db,rot) obj.OdeMomentsAdjoint(t,Y,db,rot), -obj.h, initialMoments, verbose);
            
            obj.zAdj = flip(z);
            obj.yOri = flip(y(:,12:end),1);
            obj.yAdj = flip(y(:,1:11),1);
            obj.kAdj = flip(kval);
        end
        
        % ode3 solver
        function [tout,yout,kval,Omat,Nmat] = ode3(obj, F, h, y0, verbose)
            if nargin == 3
                verbose = false;
            end
            
            if verbose
                if h > 0
                    fprintf('Integrating forward moment equations ... \n');
                else
                    fprintf('Integrating reverse moment+adjoint equations ... \n');
                end
                %msg = fprintf(['[',repelem('=',1,0),'>',repelem(' ',1,50),']']);
            end
            
            if h > 0
                intlattice = obj.lattice;
            else
                intlattice = flip(obj.lattice,1);
                intlattice(:,[1,2]) = flip(intlattice(:,[1,2]),2);
            end
            
            % variables that will be returned at each integration step
            kval = [];
            Omat = {};
            Nmat = {};
            yout = [];
            tout = [];
            
            [M,N] = size(intlattice);
            
            ii = 1;
            perc = 0; percIncrease = 2;
            % integrate for each element in the lattice
            for j = 1:M
                
                if verbose
                    progress = round( 100 * j / M );
                    if progress >= perc
                        %msg = fprintf(repmat('\b',1,msg));
                        %msg = fprintf(['[',repelem('=',1,progress),'>',repelem(' ',1,50-progress),']']);
                        fprintf(['=']);
                        perc = perc + percIncrease;
                    end
                end
                
                % start stop strength of element
                t0 = intlattice(j,1);
                t1 = intlattice(j,2);
                if ( h > 0 )
                    tsteps = linspace(t0,t1,h); %t0:h:(t1-h);
                else
                   tsteps = linspace(t1,t0,abs(h)); 
                end
                tsteps = tsteps(1:end-1);
                stepSize = abs(tsteps(2) - tsteps(1));
                db = intlattice(j,3);
                rot = intlattice(j,4);
                
                tout = [tout, tsteps]; % independent variable
                kval = [kval, zeros(1,length(tsteps))+db]; % quad strength
                
                % output variables for this element
                ytmp = zeros(length(y0),length(tsteps)+1);
                Otmp = cell(length(tsteps),1);
                Ntmp = cell(length(tsteps),1);
                
                % initial conditions are the very last set of points integrated in the previous element (except for the starting element)
                if j == 1
                    ytmp(:,1) = y0;
                else
                    ytmp(:,1) = yout(:,end);
                end
                
                % run rk3 ode solver algorithm through the element
                y = ytmp(:,1);
                for jj = 1:length(tsteps)
                    t = tsteps(jj);
                    [val,~,t2,t3] = F(t,y,db,rot);
                    s1 = stepSize.*val;
                    s2 = stepSize.*F(t+stepSize/2, y+s1./2, db, rot);
                    s3 = stepSize.*F(t+stepSize, y-s1+2*s2, db, rot);
                    y = y + (s1 + 4*s2 + s3)./6;
                    ytmp(:,jj+1) = y;
                    
                    Otmp{jj} = t2;
                    Ntmp{jj} = t3;
                    ii = ii + 1;
                end
                
                if j == 1
                    yout = [yout, ytmp];
                else
                    yout = [yout, ytmp(:,2:end)];
                end
                
                Omat = [Omat; Otmp];
                Nmat = [Nmat; Ntmp];
            end
            
            % final eval
            tout = [tout, t1];
            [val,~,t2,t3] = F(t1,yout(:,end),db,rot);
            Omat = [Omat; t2];
            Nmat = [Nmat; t3];
            kval = [kval, db];
            
            tout = tout';
            yout = yout';
            kval = kval' / obj.rigidity;
            fprintf('\n')
            
        end
        
        % moment equations being solved for a given z position and list of moments
        % Y = [Q,P,E,L]
        % Y = [Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L, phi]
        %     [   1       2    3     4       5     6    7      8     9  10 11
        function [dYdt,k_quad,O_mat,N_mat] = OdeMoments(obj, z, Y, db, rot)
            
            k_quad = db / obj.rigidity;
            psi = rot;
            k_sol = 0.0;
            
            cq = cos(2*Y(11)-2*psi);
            sq = sin(2*Y(11)-2*psi);
            
            % space charge stuff
            Q_delta = sqrt( Y(1)^2 - Y(2)^2 - Y(3)^2 );
            ab4 = 1 / Q_delta; % the 4/ab term in equation
            ca_ab4 = -Y(2) / ( (Y(1)+Q_delta)*Q_delta ); % 4c_alpha/ab
            sa_ab4 = -Y(3) / ( (Y(1)+Q_delta)*Q_delta ); % 4s_alpha/ab
            
            % Calculate O and N matrix stuff
            O_mat = [-k_sol^2/2.0 + ab4*obj.perveance, 2*k_quad*cq + ca_ab4*obj.perveance, -2*k_quad*sq + sa_ab4*obj.perveance;
                2*k_quad*cq + ca_ab4*obj.perveance, -k_sol^2/2.0 + ab4*obj.perveance, 0;
                -2*k_quad*sq + sa_ab4*obj.perveance, 0, -k_sol^2/2.0 + ab4*obj.perveance];
            
            N_mat = [0; 2*k_quad*sq - sa_ab4*obj.perveance; 2*k_quad*cq + ca_ab4*obj.perveance];
            
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
        
        function [dYdt,k_quad,O_mat,N_mat] = OdeMomentsAdjoint(obj, z, Yt, db, rot)
            Y = Yt(1:11);
            Y2 = Yt(12:end);
            
            k_quad = db / obj.rigidity;
            psi = rot;
            k_sol = 0.0;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %O & N matrix parameters calculated using original envelope equations,
            % not the adjoint version
            
            cq = cos(2*Y2(11)-2*psi);
            sq = sin(2*Y2(11)-2*psi);
            
            % space charge stuff
            
            Q_delta = sqrt( Y2(1)^2 - Y2(2)^2 - Y2(3)^2 );
            ab4 = 1 / Q_delta; % the 4/ab term in equation
            ca_ab4 = -Y2(2) / ( (Y2(1)+Q_delta)*Q_delta ); % 4c_alpha/ab
            sa_ab4 = -Y2(3) / ( (Y2(1)+Q_delta)*Q_delta ); % 4s_alpha/ab
            
            % Calculate O and N matrix stuff
            O_mat = [-k_sol^2/2.0 + ab4*obj.perveance, 2*k_quad*cq + ca_ab4*obj.perveance, -2*k_quad*sq + sa_ab4*obj.perveance;
                2*k_quad*cq + ca_ab4*obj.perveance, -k_sol^2/2.0 + ab4*obj.perveance, 0;
                -2*k_quad*sq + sa_ab4*obj.perveance, 0, -k_sol^2/2.0 + ab4*obj.perveance];
            
            N_mat = [0; 2*k_quad*sq - sa_ab4*obj.perveance; 2*k_quad*cq + ca_ab4*obj.perveance];
            
            % calculate special matrix using original position data
            [Mq,Mp,Mn] = obj.CalcAdjointMatrices(Y2);
            
            dedot1 =  Y(4:6)'*Mq + Y(10)*Y2(1:3)'*Mn - Y(1:3)'*Mp - Y(1:3)'*Mn*Y2(10);
            if 1
                % Fp 1 variable
                dedot2 = zeros(3,1);
                dqdot = zeros(3,1);
                dpdot = zeros(3,1);
                dldot = 0.0;
            else
                % Fe 2 variable
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% System of 20 equations to solve
            
            % solve regular envelope equations first
            dY1dt_1 = [ Y2(4);
                Y2(5);
                Y2(6) ];
            
            dY2dt_1 = Y2(7:9) + O_mat*Y2(1:3);
            
            dY3dt_1 = O_mat*Y2(4:6) + N_mat*Y2(10);
            
            dY4dt_1 = -N_mat'*Y2(1:3);
            
            dY5dt_1 = -1*k_sol/2.0;
            
            % solve adjoint equations
            % dQ/dz
            dY1dt_2 = Y(4:6) + dqdot;
            
            % dP/dz
            dY2dt_2 = Y(7:9) + O_mat*Y(1:3) + dpdot;
            
            % dE/dz
            dY3dt_2 = O_mat*Y(4:6) + N_mat*Y(10) + dedot1' + dedot2;
            
            % dL/dz
            dY4dt_2 = -N_mat'*Y(1:3) + dldot;
            
            % d phi / dz
            dY5dt_2 = -1*k_sol/2.0;
            
            % put all together
            dYdt2_1 = [dY1dt_1; dY2dt_1; dY3dt_1; dY4dt_1; dY5dt_1];
            dYdt2_2 = [dY1dt_2; dY2dt_2; dY3dt_2; dY4dt_2; dY5dt_2];
            dYdt = [dYdt2_2; dYdt2_1];
            
        end
        
        function [O_mat,N_mat,ACT] = CalcONandACTMatrix(obj)
            
            O_mat = cell(length(obj.z),1);
            N_mat = cell(length(obj.z),1);
            ACT = zeros(13,length(obj.z));
            
            for i = 1:length(obj.z)
                
                z = obj.z(i);
                Y = obj.y(i,:);
                Yadj = obj.yAdj(i,:);
                % find out where we are in the lattice
                % find where we are in lattice
                db = 0;
                rot = 0;
                [N,~] = size(obj.lattice);
                for j = 1:N
                    if z < obj.lattice(j,1)
                        break;
                    else
                        if z <= obj.lattice(j,2)
                            % inside quad
                            db = obj.lattice(j,3);
                            rot = obj.lattice(j,4);
                            break;
                        end
                    end
                end
                
                k_quad = db / obj.rigidity;
                psi = rot;
                k_sol = 0.0;
                
                cq = cos(2*Y(11)-2*psi);
                sq = sin(2*Y(11)-2*psi);
                
                % space charge stuff
                Q_delta = sqrt( Y(1)^2 - Y(2)^2 - Y(3)^2 );
                ab4 = 1 / Q_delta; % the 4/ab term in equation
                ca_ab4 = -Y(2) / ( (Y(1)+Q_delta)*Q_delta ); % 4c_alpha/ab
                sa_ab4 = -Y(3) / ( (Y(1)+Q_delta)*Q_delta ); % 4s_alpha/ab
                
                % Calculate O and N matrix stuff
                O_mat{i} = [-k_sol^2/2.0 + ab4*obj.perveance, 2*k_quad*cq + ca_ab4*obj.perveance, -2*k_quad*sq + sa_ab4*obj.perveance;
                    2*k_quad*cq + ca_ab4*obj.perveance, -k_sol^2/2.0 + ab4*obj.perveance, 0;
                    -2*k_quad*sq + sa_ab4*obj.perveance, 0, -k_sol^2/2.0 + ab4*obj.perveance];
                
                N_mat{i} = [0; 2*k_quad*sq - sa_ab4*obj.perveance; 2*k_quad*cq + ca_ab4*obj.perveance];
                
                % adjoint changing terms
                [Mq,Mp,Mn] = obj.CalcAdjointMatrices(Y);
                
                dedot1 =  Yadj(4:6)*Mq + Yadj(10)*Y(1:3)*Mn - Yadj(1:3)*Mp - Yadj(1:3)*Mn*Y(10);
                if 1
                    % Fp 1 variable
                    dedot2 = zeros(3,1);
                    dqdot = zeros(3,1);
                    dpdot = zeros(3,1);
                    dldot = 0.0;
                    ACT(:,i) = [dqdot;dpdot;dedot1';dedot2;dldot];
                else
                    % Fe 2 variable
                    
                end
            end
        end
        
        function [Mq,Mp,Mn] = CalcAdjointMatrices(obj, Y)
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
            Q_deltaplus = Q_delta + Y(1);
            
            V1 = -(obj.perveance/(Q_delta^2))*[Y(4)-Y(2)*Y(5)/Q_deltaplus-Y(3)*Y(6)/Q_deltaplus,
                -Y(2)*Y(4)/Q_deltaplus+Y(5),
                -Y(3)*Y(4)/Q_deltaplus+Y(6)];
            U1t = (1/Q_delta)*[Y(1), -Y(2), -Y(3)];
            
            V2 = (obj.perveance/(Q_delta*(Q_deltaplus^2)))*[Y(2)*Y(5)+Y(3)*Y(6),Y(2)*Y(4),Y(3)*Y(4)]';
            U2t = U1t + [1,0,0];
            
            V3 = -(obj.perveance/(Q_delta*Q_deltaplus))*[Y(5), Y(4), 0]';
            V4 = -(obj.perveance/(Q_delta*Q_deltaplus))*[Y(6), 0, Y(4)]';
            U3t = [0, 1, 0];
            U4t = [0,0,1];
            
            Mp = V1*U1t + V2*U2t + V3*U3t + V4*U4t;
            %%%%
            
            W1 = -(obj.perveance/Q_delta)*[1, Y(2)/Q_deltaplus, Y(3)/Q_deltaplus]';
            W2 = (obj.perveance/(Q_delta*(Q_deltaplus^2)))*[Y(2)^2+Y(3)^2,Y(2)*Y(1),Y(3)*Y(1)]';
            W3 = -(obj.perveance/(Q_delta*Q_deltaplus))*[Y(2),Y(1),0]';
            W4 = -(obj.perveance/(Q_delta*Q_deltaplus))*[Y(3),0,Y(1)]';
            
            Mq = W1*U1t + W2*U2t + W3*U3t + W4*U4t;
            %%%%
            
            X1 = -obj.perveance*[0, Y(3)/(Q_deltaplus*(Q_delta^2)), -Y(2)/(Q_deltaplus*(Q_delta^2))]';
            X2 = -obj.perveance*[0, Y(3)/(Q_delta*(Q_deltaplus^2)), -Y(2)/(Q_delta*(Q_deltaplus^2))]';
            X3 = obj.perveance*[0,0,-1/(Q_delta*Q_deltaplus)]';
            X4 = obj.perveance*[0,1/(Q_delta*Q_deltaplus),0]';
            
            Mn = X1*U1t + X2*U2t + X3*U3t + X4*U4t;
            
        end
        
        function [intVal] = AdjointIntegral(obj,O_opt,N_opt)
            %FIND_INT Summary of this function goes here
            
            % calculate adjoint integral broken up into 4 pieces
            int1 = zeros(1,length(obj.z));
            int2 = zeros(1,length(obj.z));
            int3 = zeros(1,length(obj.z));
            int4 = zeros(1,length(obj.z));
            
            for i = 1:length(obj.z)
                int1(i) = obj.yAdj(i,4:6)*(O_opt{i}*obj.y(i,1:3)');
                int2(i) = obj.yAdj(i,10)*(obj.y(i,1:3)*N_opt{i});
                int3(i) = -1*obj.yAdj(i,1:3)*(O_opt{i}*obj.y(i,4:6)');
                int4(i) = -1*obj.yAdj(i,1:3)*N_opt{i}*obj.y(i,10);
            end
            
            intVal = trapz(obj.z,int1+int2+int3+int4);
        end
        
        function [FoM, FoMp, dQy, dPy, dEy, dLy] = GetFAndDF1(obj)
            
            y = obj.y;
            k0 = obj.k0;
            
            % match location
            mloc = (length(obj.z));
            mloc2 = (length(obj.z)-1)/2.0;
            
            % calculate adjoint related variables
            dPy = zeros(3,1); gPy = dPy;
            dEy = zeros(3,1); gEy = dEy;
            dQy = zeros(3,1); gQy = dQy;
            dLy = zeros(1,1); gLy = dLy;
            
            % figure of merit periodic
            FoM1 = 0.5 * k0^2 * sum((y(mloc,1:3) - y(1,1:3)).^2);
            FoM2 = 0.5 * sum((y(mloc,4:6) - y(1,4:6)).^2);
            FoM3 = 0.5 * k0^(-2) * sum((y(mloc,7:9) - y(1,7:9)).^2);
            FoM4 = 0.5 * sum((y(mloc,10) - y(1,10)).^2);
            
            % dF/dQ
            gQy(1:3) = k0^2 * (y(mloc,1:3) - y(1,1:3));
            % dF/dP
            gPy(1:3) = (y(mloc,4:6) - y(1,4:6));
            % dF/dE
            gEy(1:3) = k0^(-2) * (y(mloc,7:9) - y(1,7:9));
            % dF/dL
            gLy(1)   = (y(mloc,10) - y(1,10));
            
            dQy = - gEy;
            dPy =   gPy;
            dEy = - gQy;
            dLy  = - gLy;
            
            FoM = (abs(FoM1) + abs(FoM2) + abs(FoM3) + abs(FoM4));
            FoMp = [abs(FoM1),abs(FoM2),abs(FoM3),abs(FoM4)];
            
        end
        
        function [df] = CalcFoMGradientX(obj)
            
            if 1
                df = zeros(11,1);
                k0 = obj.k0;
                
                % Q
                df(1:3) = -k0^(-1) * (obj.yAdj(1,7:9) - obj.yAdj(end,7:9));
                % P
                df(4:6) = (obj.yAdj(1,4:6) - obj.yAdj(end,4:6));
                % E
                df(7:9) = -k0 * (obj.yAdj(1,1:3) - obj.yAdj(end,1:3));
                % L
                df(10) = -(obj.yAdj(1,10) - obj.yAdj(end,10));
                % phi
                df(11) = 0.0; % solenoid phi larmor angle
            else
                df = zeros(9,1);
                k0 = obj.k0;
                
                % Q
                df(1:3) = -k0^(-1) * (obj.yAdj(1,7:9) - obj.yAdj(end,7:9));
                % P
                df(4:6) = (obj.yAdj(1,4:6) - obj.yAdj(end,4:6));
                % E
                df(7) = -k0 * (obj.yAdj(1,3) - obj.yAdj(end,3));
                % L
                df(8) = -(obj.yAdj(1,10) - obj.yAdj(end,10));
                % phi
                df(9) = 0.0; % solenoid phi larmor angle
            end
            
        end
        
        function [djs] = CalcFoMGradientJ(obj)
            idx = 1;
            % j terms term
            % ex2
            j1 = (obj.y(idx,7)+obj.y(idx,8)) .* (obj.y(idx,1)+obj.y(idx,2)) - 0.5*(obj.y(idx,4)+obj.y(idx,5)).^2;
            dj1 = zeros(11,1);         
            dj1(1) = (obj.y(idx,7)+obj.y(idx,8));
            dj1(2) = (obj.y(idx,7)+obj.y(idx,8));
            dj1(3) = 0.0;
            dj1(4) = -(obj.y(idx,4)+obj.y(idx,5));
            dj1(5) = -(obj.y(idx,4)+obj.y(idx,5));
            dj1(6) = 0.0;
            dj1(7) = (obj.y(idx,1)+obj.y(idx,2));
            dj1(8) = (obj.y(idx,1)+obj.y(idx,2));
            dj1(9) = 0.0;
            dj1(10) = 0.0;
            dj1(11) = 0.0;
            % ey2
            j2 = (obj.y(idx,7)-obj.y(idx,8)) .* (obj.y(idx,1)-obj.y(idx,2)) - 0.5*(obj.y(idx,4)-obj.y(idx,5)).^2;
            dj2 = zeros(11,1);             
            dj2(1) = (obj.y(idx,7)-obj.y(idx,8));
            dj2(2) = -(obj.y(idx,7)-obj.y(idx,8));
            dj2(3) = 0.0;
            dj2(4) = -(obj.y(idx,4)-obj.y(idx,5));
            dj2(5) = (obj.y(idx,4)-obj.y(idx,5));
            dj2(6) = 0.0;
            dj2(7) = (obj.y(idx,1)-obj.y(idx,2));
            dj2(8) = -(obj.y(idx,1)-obj.y(idx,2));
            dj2(9) = 0.0;
            dj2(10) = 0.0;
            dj2(11) = 0.0;            
            % com
            dj3 = zeros(11,1);
            % Q
            dj3(1) = obj.y(idx,7);
            dj3(2) = obj.y(idx,8);
            dj3(3) = obj.y(idx,9);
            % P
            dj3(4) = -obj.y(idx,4);
            dj3(5) = -obj.y(idx,5);
            dj3(6) = -obj.y(idx,6);
            % E
            dj3(7) = obj.y(idx,1);
            dj3(8) = obj.y(idx,2);
            dj3(9) = obj.y(idx,3);
            % L
            dj3(10) = obj.y(idx,10);
            % phi
            dj3(11) = 0; % solenoid phi larmor angle                
            
            djs = zeros(11,2);
            djs(:,1) = dj1;
            djs(:,2) = dj2;
            %djs(:,3) = dj3;
            
        end        
        
        function PlotBeamSize(obj)
            
            figure; hold on;
            plot(obj.z,obj.y(:,1)+obj.y(:,2));
            plot(obj.z,obj.y(:,1)-obj.y(:,2));
            
            kquad = obj.k / max(obj.k);
            kquad = kquad * max(obj.y(:,1)) * 0.10;
            plot(obj.z,kquad,'k-','Linewidth',1);
            
            xlabel('Z position (m)'); ylabel('Moments [m^2]'); title('Beam size through lattice');
            legend('\langle x^2 \rangle','\langle y^2 \rangle','K_{quad}','Location','northeastoutside');
            grid on;
        end
        
        function PlotBeamEmit(obj)
            ex = ( (obj.y(:,7)+obj.y(:,8)) .* (obj.y(:,1)+obj.y(:,2)) - 0.5*(obj.y(:,4)+obj.y(:,5)).^2 );
            ey = ( (obj.y(:,7)-obj.y(:,8)) .* (obj.y(:,1)-obj.y(:,2)) - 0.5*(obj.y(:,4)-obj.y(:,5)).^2 );
            
            figure; hold on;
            plot(obj.z,ex);
            plot(obj.z,ey);
            
            kquad = obj.k / max(obj.k);
            kquad = kquad * max(ex) * 0.10;
            plot(obj.z,kquad,'k-','Linewidth',1);
            
            xlabel('Z position (m)'); ylabel('Emittance [m*rad]'); title('Beam emittance through lattice');
            legend('\epsilon_x^2','\epsilon_y^2','K_{quad}','Location','northeastoutside');
            grid on;            
        end
        
    end
end

