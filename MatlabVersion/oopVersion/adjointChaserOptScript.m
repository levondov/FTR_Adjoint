endOpt = false;

runXopt = true;
runYopt = true;

% generate some random numbers
rng('default');
rng(2); % use this to repeat runs, seed = 2 etc.

% Initial conditions
% 0mA
X0 = [
   0.000407212723699
  -0.000183484865165
                   0
   0.000000005174747
   0.000000006065963
                   0
   0.118445742459199
   0.053370142482708
                   0
                   0
                   0] * 1e-3;  

for i = 1:length(X0)
    X0(i) = X0(i) + X0(i)*rand*0.05;
end
               
momOpt = adjointChaserOpt(0, 2, X0, X0);
momOpt.useJConstraint = true;
momOpt.lambdaY = 1.0e-12;
momOpt.momY.M = 2;
momOpt.momY.QTE = 3.933066132264529e-07;

momOpt = momOpt.initX();
momOpt = momOpt.initY();

momOpt.momY.PlotBeamSize

figure;
plot(momOpt.momY.z, momOpt.momY.y(:,1),'linewidth',2); hold on;
plot(momOpt.momY.z, repmat(momOpt.momY.QTE,length(momOpt.momY.z),1),'linewidth',2);
xlabel('Z [m]');
ylabel('Moments');
legend('Q+','QTE');
title('Final Solution, Q+ vs z with QTE');

if (runYopt) 
    momOpt = momOpt.findStartingGammaY(); 
end
if (runXopt) 
    momOpt = momOpt.findStartingGammaX(); 
end

%% Gradient descent algorithm
% Here we run our gradient descent algorithm
Xntmp = [];
while ~endOpt % let it run forever and break with ctrl-c , or if conditions are met at the end of this while loop it will break automatically
    ii = 0; iiy = 0; iix = 0;

    if ( runYopt && runXopt )
        % if this is false, we stop taking steps and recompute gradients
        yStepCondition = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) < (momOpt.WY_h(end-1) + momOpt.lambdaY*momOpt.FeY_h(end-1));
        xStepCondition = (momOpt.WX_h(end)) < (momOpt.WX_h(end-1));        
        ii = 1;
        while ( true == yStepCondition && true == xStepCondition )
            fprintf(['Iterating ',num2str(iiy),'\n']);

            % take steps
            momOpt = momOpt.takeStepY();
            momOpt = momOpt.takeStepX();
            ii = ii + 1;

            % iterative scheme to get J back to starting J0 value
            momOpt = momOpt.iterateJForY();  
            momOpt = momOpt.iterateJForX(); 

            % Run moment equations and update FoM
            momOpt = momOpt.UpdateFoMAll();

            % if we take 50+ steps without having to recompute gradients, our
            % steps are probably too small, let's start increasing gamma
            if ( ii > 50 )
                % increase gamma
                fprintf(['Too many small steps, increasing gamma... \n']);
                momOpt.gammaY_h(end+1) = momOpt.gammaY_h(end) * 2.0;
                momOpt.gammaX_h(end+1) = momOpt.gammaX_h(end) * 2.0;
                ii = 0;
            end

            % update best values
            if (momOpt.WY_h(end) < momOpt.WbestY_h(end))
               momOpt.WbestY_h(end+1) = momOpt.WY_h(end); 
            end
            if (momOpt.FeY_h(end) < momOpt.FebestY_h(end))
               momOpt.FebestY_h(end+1) = momOpt.FeY_h(end); 
            end  
            if ( (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) <  momOpt.WandFebestY_h(end) )
                momOpt.WandFebestY_h(end+1) = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end));
            end            
            if (momOpt.WX_h(end) < momOpt.WbestX_h(end))
               momOpt.WbestX_h(end+1) = momOpt.WX_h(end); 
            end

            % if this is false, we stop taking steps and recompute gradients
            yStepCondition = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) < (momOpt.WY_h(end-1) + momOpt.lambdaY*momOpt.FeY_h(end-1));
            xStepCondition = (momOpt.WX_h(end)) < (momOpt.WX_h(end-1));          

        end
    elseif ( true == runYopt )
         % if this is false, we stop taking steps and recompute gradients
        yStepCondition = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) < (momOpt.WY_h(end-1) + momOpt.lambdaY*momOpt.FeY_h(end-1));                   
        iiy=1;
        % while the FoM keeps decreasing
        while ( true == yStepCondition )
            fprintf(['Iterating ',num2str(iiy),'\n']);
    
            momOpt = momOpt.takeStepY();
    
            iiy = iiy + 1;
    
            % iterative scheme to get J back to starting J0 value
            momOpt = momOpt.iterateJForY();
            
            % Run moment equations
            momOpt = momOpt.UpdateFoMY();
            
            % if we take 50+ steps without having to recompute gradients, our
            % steps are probably too small, let's start increasing gamma
            if ( iiy > 50 )
                % increase gamma
                fprintf(['Too many small steps, increasing gamma... \n']);
                momOpt.gammaY_h(end+1) = momOpt.gammaY_h(end) * 2.0;
                iiy = 0;
            end
            
            if (momOpt.WY_h(end) < momOpt.WbestY_h(end))
               momOpt.WbestY_h(end+1) = momOpt.WY_h(end); 
            end
            if (momOpt.FeY_h(end) < momOpt.FebestY_h(end))
               momOpt.FebestY_h(end+1) = momOpt.FeY_h(end); 
            end  
            if ( (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) <  momOpt.WandFebestY_h(end) )
                momOpt.WandFebestY_h(end+1) = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end));
            end

            yStepCondition = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) < (momOpt.WY_h(end-1) + momOpt.lambdaY*momOpt.FeY_h(end-1));                   
                    
        end
    elseif ( true == runXopt )
        xStepCondition = (momOpt.WX_h(end)) < (momOpt.WX_h(end-1));
        iix=1;
        % while the FoM keeps decreasing
        while ( true == xStepCondition )
            fprintf(['Iterating ',num2str(iix),'\n']);
    
            momOpt = momOpt.takeStepX();
    
            iix = iix + 1;
    
            % iterative scheme to get J back to starting J0 value
            momOpt = momOpt.iterateJForX();
            
            % Run moment equations
            momOpt = momOpt.UpdateFoMX();
            
            % if we take 50+ steps without having to recompute gradients, our
            % steps are probably too small, let's start increasing gamma
            if ( iix > 50 )
                % increase gamma
                fprintf(['Too many small steps, increasing gamma... \n']);
                momOpt.gammaX_h(end+1) = momOpt.gammaX_h(end) * 2.0;
                iix = 0;
            end
            
            if (momOpt.WX_h(end) < momOpt.WbestX_h(end))
               momOpt.WbestX_h(end+1) = momOpt.WX_h(end); 
            end

            xStepCondition = (momOpt.WX_h(end)) < (momOpt.WX_h(end-1));               
        end
    else
        % not running any optimization??
    end 

    % if FoM is no longer decreasing, lets recalculate the adjoint
    % equations
    
    % recompute adjoint equation stuff for new direction
    fprintf(['Recomputing adjoint equations \n']);
    
    % grab last good settings
    momOpt = momOpt.GrabLastGoodValues();
    
    % calc adjoint equations
    momOpt = momOpt.RecalculateGradients();

    % this should be called if we are running both X and Y opt
    if (ii == 2) % meaning no improving from recalculating gradient.
    % recalculating gives no improvement, lets try adjusting gamma
        %momOpt = momOpt.findStartingGammaY();
        %momOpt = momOpt.findStartingGammaX();        
        
        if (momOpt.WY_h(end) < momOpt.WbestY_h(end))
           momOpt.WbestY_h(end+1) = momOpt.WY_h(end); 
        end  
        if (momOpt.FeY_h(end) < momOpt.FebestY_h(end))
           momOpt.FebestY_h(end+1) = momOpt.FeY_h(end); 
        end   
        if ( (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) <  momOpt.WandFebestY_h(end) )
            momOpt.WandFebestY_h(end+1) = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end));
        end    
        if (momOpt.WX_h(end) < momOpt.WbestX_h(end))
           momOpt.WbestX_h(end+1) = momOpt.WX_h(end); 
        end            
    end
    
    % this should only be called if we are running just Y opt
    if (iiy == 2) % meaning no improving from recalculating gradient.
        % recalculating gives no improvement, lets try adjusting gamma
        %momOpt = momOpt.findStartingGammaY();
        
        if (momOpt.WY_h(end) < momOpt.WbestY_h(end))
           momOpt.WbestY_h(end+1) = momOpt.WY_h(end); 
        end  
        if (momOpt.FeY_h(end) < momOpt.FebestY_h(end))
           momOpt.FebestY_h(end+1) = momOpt.FeY_h(end); 
        end   
        if ( (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end)) <  momOpt.WandFebestY_h(end) )
            momOpt.WandFebestY_h(end+1) = (momOpt.WY_h(end) + momOpt.lambdaY*momOpt.FeY_h(end));
        end        
    end

    % this should only be called ifw e are running just X opt
    if (iix == 2) % meaning no improving from recalculating gradient.
        % recalculating gives no improvement, lets try adjusting gamma
        %momOpt = momOpt.findStartingGammaX();
        
        if (momOpt.WX_h(end) < momOpt.WbestX_h(end))
           momOpt.WbestX_h(end+1) = momOpt.WX_h(end); 
        end        
    end    
    
    if momOpt.WX_h(end) < 1e-30 % break if we get this low
        endOpt = true;
    end
    if length(momOpt.WX_h(end)) > 100000 % break if we have done this many iterations
        endOpt = true;
    end
end

%% Some plots to post afterwards;

momOpt.momY.PlotBeamSize
momOpt.momX.PlotBeamSize

figure;
subplot(1,3,1); hold on;
plot(log10(momOpt.WbestY_h(2:end)),'linewidth',2); title('W(Y,a) & W(X,a)'); ylabel('log10(W)');
plot(log10(momOpt.WbestX_h(2:end)),'linewidth',2); legend('W(Y,a)','W(X,a)')
subplot(1,3,2); plot(log10(momOpt.FebestY_h),'linewidth',2); title('FE (Y,a)'); ylabel('log10(FE)'); xlabel('Iterations');
subplot(1,3,3); hold on;
plot(log10(momOpt.WandFebestY_h),'linewidth',2); title('W + \lambda FE'); ylabel('log10(W+FE)');
plot(log10(momOpt.WbestX_h(2:end)),'linewidth',2); legend('W(Y,a)','W(X,a)')

figure; 
subplot(1,3,1); hold on; 
plot(log10(momOpt.WY_h),'linewidth',2); ylabel('log10(W)'); xlabel('Iterations'); title('W vs iterations');
plot(log10(momOpt.WX_h),'linewidth',2);
legend('W(Y,a)','W(X,a)');
subplot(1,3,2); plot(log10(momOpt.FeY_h)); ylabel('log10(FE)'); xlabel('Iterations'); title('FE vs iterations');
subplot(1,3,3); plot(log10(momOpt.WY_h + momOpt.lambdaY * momOpt.FeY_h)); ylabel('log10(FE)'); xlabel('Iterations'); title('W + \lambda FE vs iterations');

figure;
plot(momOpt.momY.z, momOpt.momY.y(:,1),'linewidth',2); hold on;
plot(momOpt.momY.z, repmat(momOpt.momY.QTE,length(momOpt.momY.z),1),'linewidth',2);
xlabel('Z [m]');
ylabel('Moments');
legend('Q+','QTE');
title('Final Solution, Q+ vs z with QTE');

figure;
subplot(1,3,1); hold on;
plot(momOpt.momY.z, abs(momOpt.momY.y(:,1) - momOpt.momX.y(:,1)) ./ momOpt.momX.y(:,1),'linewidth',2);
xlabel('Z [m]')
ylabel('Moments [m^2]');
title('Q+ comparison');
subplot(1,3,2); hold on;
plot(momOpt.momY.z, abs(momOpt.momY.y(:,4) - momOpt.momX.y(:,4)) ./ momOpt.momX.y(:,4),'linewidth',2);
xlabel('Z [m]')
ylabel('Moments [m^2]');
title('P+ comparison');
subplot(1,3,3); hold on;
plot(momOpt.momY.z, abs(momOpt.momY.y(:,7) - momOpt.momX.y(:,7)) ./ momOpt.momX.y(:,7),'linewidth',2);
xlabel('Z [m]')
ylabel('Moments [m^2]');
title('E+ comparison');

%%

figure; 
xx = momOpt.Xn_h(:,1:6918);
yy = -momOpt.Yn_h(:,1:end);
diffx = xx-yy;

plot(vecnorm(diffx),'linewidth',2)

