classdef adjointChaserOpt
    %ADJOINTCHASEROPT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        momX;
        momY;

        % X related data
        WX_h;
        WpX_h;
        Xn_h;
        gammaX_h;
        dWX_h;
        djX_h;
        jX_h;
        WbestX_h;

        % Y related data
        WY_h;
        WpY_h;
        Yn_h;
        gammaY_h;
        dWY_h;
        djY_h;
        jY_h;
        WbestY_h;
        FeY_h;
        FebestY_h;
        dFeY_h;
        lambdaY = 1.0;
        WandFebestY_h;

        % opt properties
        useJConstraint = true; % include the J constraint in optimization
    end

    methods
        function obj = adjointChaserOpt(current, period, X0, Y0)
            % calculate E+ and E- for a give emittance
            ex2 = 6.2e-6^2; % target emittance X
            ey2 = 6.2e-6^2; % target emittance Y
            denom1 = (X0(1)+X0(2)); denom2 = (X0(1)-X0(2));
            numer1 = ex2 + 0.5*(X0(4)+X0(5)).^2; numer2 = ey2 + 0.5*(X0(4)-X0(5)).^2;
            Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
            Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
            X0(7) = Ep;
            X0(8) = Em;
            denom1 = (Y0(1)+Y0(2)); denom2 = (Y0(1)-Y0(2));
            numer1 = ex2 + 0.5*(Y0(4)+Y0(5)).^2; numer2 = ey2 + 0.5*(Y0(4)-Y0(5)).^2;
            Ep = 0.5*((numer1 / denom1) + (numer2 / denom2));
            Em = 0.5*((numer1 / denom1) - (numer2 / denom2));
            Y0(7) = Ep;
            Y0(8) = Em;
            Xn = X0';
            YN = X0';

            % setup moment object
            momX = MomentSolverPeriodicChaser(10e3, current, X0); % energy, beam current, initial conditions
            momX.h = 100; % integration steps in each element
            momY = MomentSolverPeriodicChaser(10e3, current, Y0); % energy, beam current, initial conditions
            momY.h = 100; % integration steps in each element

            % here we have 32 quads each with a adjustable strength dB, thus our a vector is length 32
            qlength = 0.02;
            qstart = [0.0, 0.03, 0.07];
            db = [0.5, -0.5, 0.5];
            qend = [qstart(1)+qlength/2.0, qstart(2)+qlength, qstart(3)+qlength/2.0];
            qrot = [0.0, 0.0, 0.0, 0.0];
            zstart = 0.0;
            zend = 0.08;

            % setup lattice
            momX = momX.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);
            momY = momY.CreateLatticeProfile(db,qstart,qend,qrot, zstart, zend, period, false);

            obj.momX = momX;
            obj.momY = momY;
        end

        function obj = initX(obj)
            % run moment + adjoint equations
            obj.momX = obj.momX.RunMoments();
            obj.momX = obj.momX.RunMomentsAdjoint();

            % Grab FoM
            [w0,w0p] = obj.momX.GetW();

            % gradient
            df0 = obj.momX.GetWGradientX();

            % grad J gradient
            [dj0,j0] = obj.momX.GetJGradientX();
            obj.jX_h = j0;
            obj.djX_h = {dj0};

            % gradient descent parameter
            gammaX = (w0/sum(df0.^2));

            % init arrays
            obj.WX_h = w0;
            obj.WpX_h = w0p;
            obj.Xn_h = obj.momX.initialMoments;
            obj.gammaX_h = gammaX;
            obj.dWX_h = df0;
            obj.WbestX_h = w0;
        end

        function obj = initY(obj)
            % run moment + adjoint equations
            obj.momY = obj.momY.RunMoments();
            obj.momY = obj.momY.RunMomentsAdjoint();

            % Grab FoM
            [w0,w0p] = obj.momY.GetW();

            % gradient
            df0 = obj.momY.GetWGradientX();

            % grad J gradient
            [dj0,j0] = obj.momY.GetJGradientX();
            obj.jY_h = j0;
            obj.djY_h = {dj0};

            % gradient descent parameter
            gammaY = (w0/sum(df0.^2));

            % init arrays
            obj.WY_h = w0;
            obj.WpY_h = w0p;
            obj.Yn_h = obj.momY.initialMoments;
            obj.gammaY_h = gammaY;
            obj.dWY_h = df0;
            obj.WbestY_h = w0;

            % Fe stuff
            obj.momY = obj.momY.RunMomentsAdjointFE();
            % Grab FoM
            [fe0] = obj.momY.GetFE();
            % gradient
            df0 = obj.momY.GetFEGradientX();

            obj.FeY_h = fe0;
            obj.FebestY_h = fe0;
            obj.dFeY_h = df0;

            obj.WandFebestY_h = w0 + obj.lambdaY * fe0;
        end

        function obj = takeStepX(obj, X0)
            % grab "a" piece
            a = obj.calcAvector(obj.djX_h{end}, obj.dWX_h(:,end), obj.gammaX_h(end));

            % sum J component
            Jpart = 0;
            djX = obj.djX_h{end};
            for i = 1:length(a)
                Jpart = Jpart + a(i)*djX(:,i);
            end

            if (obj.useJConstraint == false)
                Jpart = 0;
            end

            if (nargin == 1)
                X0 =  obj.Xn_h(:,end);
            end
            % take step
            X = X0 - obj.gammaX_h(end) * obj.dWX_h(:,end) ... % dW
                - Jpart;

            obj.Xn_h(:,end+1) = X;
        end

        function obj = takeStepY(obj, Y0)

            % grab "a" piece for gradient of J compoenent
            a = obj.calcAvector(obj.djY_h{end}, obj.dWY_h(:,end), obj.gammaY_h(end));

            % sum J component
            Jpart = 0;
            djY = obj.djY_h{end};
            for i = 1:length(a)
                Jpart = Jpart + a(i)*djY(:,i);
            end

            if (obj.useJConstraint == false)
                Jpart = 0;
            end

            if (nargin == 1)
                Y0 =  obj.Yn_h(:,end);
            end

            % take step
            Y = Y0 - obj.gammaY_h(end) * obj.dWY_h(:,end) ... % dW
                - obj.gammaY_h(end) * obj.lambdaY * obj.dFeY_h(:,end)... % dFe
                - Jpart; % dJ

            obj.Yn_h(:,end+1) = Y;
        end

        function obj = findStartingGammaX(obj)
            % repeat this step until we get a FoM lower than what we started with
            fprintf(['FoM: ',num2str(obj.WX_h(end)),'\n']);
            ii = 0;
            f0 = obj.WX_h(end);
            X0 = obj.Xn_h(:,end);
            while obj.WX_h(end) >= f0
                obj.gammaX_h(end+1) = obj.gammaX_h(end)/2.0; % keep updating gamma here

                % Take a step
                obj = obj.takeStepX(X0);

                % Run moment equations
                obj.momX.initialMoments = obj.Xn_h(:,end)';
                obj.momX = obj.momX.RunMoments();

                % calculate FoM
                [obj.WX_h(end+1), obj.WpX_h(end+1,:)] = obj.momX.GetW();
                fprintf(['FoM: ',num2str(obj.WX_h(end)),'\n']);
                ii = ii +1;
                if ( ii > 40 )
                    fprintf(['Can not take a first step to lower FoM W \n']);
                    break;
                end
            end
            obj.WbestX_h(end+1) = obj.WX_h(end);
        end

        function obj = findStartingGammaY(obj)
            % repeat this step until we get a FoM lower than what we started with
            fprintf(['FoM: ',num2str(obj.WY_h(end) + obj.lambdaY * obj.FeY_h(end)),'\n']);
            ii = 0;
            f0 = obj.WY_h(end) + obj.lambdaY * obj.FeY_h(end);
            Y0 = obj.Yn_h(:,end);
            while (obj.WY_h(end) + obj.lambdaY * obj.FeY_h(end)) >= f0
                obj.gammaY_h(end+1) = obj.gammaY_h(end)/2.0; % use same thing for Y

                % Take a step
                obj = obj.takeStepY(Y0);

                % Run moment equations
                obj.momY.initialMoments = obj.Yn_h(:,end)';
                obj.momY = obj.momY.RunMoments();

                % calculate FoM
                [obj.WY_h(end+1), obj.WpY_h(end+1,:)] = obj.momY.GetW();
                [obj.FeY_h(end+1)] = obj.momY.GetFE();
                fprintf(['FoM: ',num2str(obj.WY_h(end)+obj.lambdaY *obj.FeY_h(end)),'\n']);
                ii = ii +1;
                if ( ii > 40 )
                    fprintf(['Can not take a first step to lower FoM W \n']);
                    break;
                end
            end
            obj.WbestY_h(end+1) = obj.WY_h(end);
        end

        function obj = UpdateFoMAll(obj)
            % update W(X,a) FoM
            obj.momX.initialMoments = obj.Xn_h(:,end)';
            obj.momX = obj.momX.RunMoments();
            [obj.WX_h(end+1), obj.WpX_h(end+1,:)] = obj.momX.GetW();

            % update W(Y,a) and FE(Y,a) FoM
            obj.momY.initialMoments = obj.Yn_h(:,end)';
            obj.momY = obj.momY.RunMoments();
            [obj.WY_h(end+1), obj.WpY_h(end+1,:)] = obj.momY.GetW();
            [obj.FeY_h(end+1)] = obj.momY.GetFE();
            fprintf(['FoM X: ',num2str(obj.WX_h(end)),' | FoM Y: ',num2str(obj.WY_h(end)),' | FoM FE Y: ',num2str(obj.FeY_h(end)),'\n']);
        end

        function obj = UpdateFoMY(obj)
            % update W(Y,a) and FE(Y,a) FoM
            obj.momY.initialMoments = obj.Yn_h(:,end)';
            obj.momY = obj.momY.RunMoments();
            [obj.WY_h(end+1), obj.WpY_h(end+1,:)] = obj.momY.GetW();
            [obj.FeY_h(end+1)] = obj.momY.GetFE();
            fprintf(['FoM Y: ',num2str(obj.WY_h(end)),' | FoM FE Y: ',num2str(obj.FeY_h(end)),'\n']);
        end

        function obj = UpdateFoMX(obj)
            % update W(X,a) FoM
            obj.momX.initialMoments = obj.Xn_h(:,end)';
            obj.momX = obj.momX.RunMoments();
            [obj.WX_h(end+1), obj.WpX_h(end+1,:)] = obj.momX.GetW();
            fprintf(['FoM X: ',num2str(obj.WX_h(end)),'\n']);
        end

        function obj = RecalculateGradients(obj)
            % X
            obj.momX.initialMoments = obj.Xn_h(:,end);
            obj.momX = obj.momX.RunMoments();
            obj.momX = obj.momX.RunMomentsAdjoint();
            % gradient
            obj.dWX_h(:,end+1) = obj.momX.GetWGradientX();

            % Y
            obj.momY.initialMoments = obj.Yn_h(:,end);
            obj.momY = obj.momY.RunMoments();
            obj.momY = obj.momY.RunMomentsAdjoint();
            obj.dWY_h(:,end+1) = obj.momY.GetWGradientX();

            [dj,j0] = obj.momY.GetJGradientX();
            obj.djY_h{end+1} = dj;
            obj.jY_h(end+1,:) = j0;           

            obj.momY = obj.momY.RunMomentsAdjointFE();
            obj.dFeY_h(:,end+1) = obj.momY.GetFEGradientX();
        end

        function obj = GrabLastGoodValues(obj)
            % grab last good settings
            obj.Xn_h(:,end+1) = obj.Xn_h(:,end-1);
            obj.Yn_h(:,end+1) = obj.Yn_h(:,end-1);   
            obj.WX_h(end+1) = obj.WX_h(end-1);
            obj.WY_h(end+1) = obj.WY_h(end-1); 
            obj.WpX_h(end+1,:) = obj.WpX_h(end-1,:);
            obj.WpY_h(end+1,:) = obj.WpY_h(end-1,:);    
            obj.FeY_h(end+1) = obj.FeY_h(end-1);
        end
        
        % Related to J calculation
        function [a] = calcAvector(obj, djs, dw, b)
            [~,N] = size(djs);
            a = zeros(N,1);

            U = obj.calcUvector(djs, dw);
            M = obj.calcMvector(djs, dw);

            % calculate a
            for i = 1:N
                sumComp = 0;
                for j = 1:N
                    sumComp = sumComp + M(i,j)*U(j);
                end
                a(i) = -b * sumComp;
            end
        end

        % Related to J calculation
        function [U] = calcUvector(obj, djs, dw)
            [~,N] = size(djs);

            U = zeros(N,1);

            for i = 1:N
                U(i) = dot( dw, djs(:,i) );
            end
        end

        % Related to J calculation
        function [M] = calcMvector(obj, djs, dw)
            [~,N] = size(djs);

            % calculate M matrix
            M1 = zeros(N,N);
            for i = 1:N
                for j = 1:N
                    M1(i,j) = dot( djs(:,i), djs(:,j) );
                end
            end
            % inverse it
            M = inv(M1);
        end

        %
        % Iterative scheme to fix J
        function obj = iterateJForY(obj)
            %
            if (obj.useJConstraint == false)
                return;
            end

            % grab starting j0, the ideal j0 that we want to have.
            j0 = obj.jY_h(1,:);
            mu = 1;
            [~,N] = size(j0);

            % do 10 iterations
            for ii = 1:3
                % calculate J and dJ/dX for the iteration
                obj.momY.initialMoments = obj.Yn_h(:,end); % grab last set of moment values
                obj.momY = obj.momY.RunMoments(); % solve moment equations
                [dj,js] = obj.momY.GetJGradientX(); % grab the resulting J and dJ/dW values for the given X moments

                M = obj.calcMvector(dj, obj.dWY_h(:,end)); % calculate M matrix
                obj.djY_h{end+1} = dj; % update history
                obj.jY_h(end+1,:) = js; % update history

                % calculate sums for the iteration scheme
                % Xn = Xn-1 + \sum a_i,n-1 * dJ_idX (see Tom's equations)
                fullVal = 0;
                for i = 1:N
                    aVal = 0;
                    for j = 1:N
                        aVal = aVal + M(i,j)*( j0(j)-js(j) );
                    end
                    fullVal = fullVal + aVal*dj(:,i);
                end

                % iterate
                obj.Yn_h(:,end+1) = obj.Yn_h(:,end) + mu*fullVal;
            end
        end
        
        function obj = iterateJForX(obj)
            %
            if (obj.useJConstraint == false)
                return;
            end

            % grab starting j0, the ideal j0 that we want to have.
            j0 = obj.jX_h(1,:);
            mu = 1;
            [~,N] = size(j0);

            % do 10 iterations
            for ii = 1:3
                % calculate J and dJ/dX for the iteration
                obj.momX.initialMoments = obj.Xn_h(:,end); % grab last set of moment values
                obj.momX = obj.momX.RunMoments(); % solve moment equations
                [dj,js] = obj.momX.GetJGradientX(); % grab the resulting J and dJ/dW values for the given X moments

                M = obj.calcMvector(dj, obj.dWX_h(:,end)); % calculate M matrix
                obj.djX_h{end+1} = dj; % update history
                obj.jX_h(end+1,:) = js; % update history

                % calculate sums for the iteration scheme
                % Xn = Xn-1 + \sum a_i,n-1 * dJ_idX (see Tom's equations)
                fullVal = 0;
                for i = 1:N
                    aVal = 0;
                    for j = 1:N
                        aVal = aVal + M(i,j)*( j0(j)-js(j) );
                    end
                    fullVal = fullVal + aVal*dj(:,i);
                end

                % iterate
                obj.Xn_h(:,end+1) = obj.Xn_h(:,end) + mu*fullVal;
            end
        end

    end % methods
end % class

