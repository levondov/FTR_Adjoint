function [z,y] = diff_in_parts(params,init_cond,reverse_flag)
%DIFF_IN_PARTS Summary of this function goes here
% integrate in parts through each section

h_drift = 0.000001;
h_quad = 0.0001;
h_sol = 0.0001;
solenoid_flag = 0;

if ~reverse_flag
    
    % Drift distance up to first quad.
    z_interval = [0.0,params(4)];
    h = h_drift;
    z1 = z_interval(1):h:z_interval(2);
    [y1] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), init_cond);
    
    % through quad 1
    z_interval = [params(4),params(4)+params(5)];
    h = h_quad;
    z2 = z_interval(1):h:z_interval(2);
    [y2] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y1(end,:)');
    
    % through drift 2
    z_interval = [params(4)+params(5),params(7)];
    h = h_drift;
    z3 = z_interval(1):h:z_interval(2);
    [y3] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y2(end,:)');
    
    % through quad 2
    z_interval = [params(7),params(7)+params(8)];
    h = h_quad;
    z4 = z_interval(1):h:z_interval(2);
    [y4] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y3(end,:)');
    
    % through drift 3
    z_interval = [params(7)+params(8),params(10)];
    h = h_drift;
    z5 = z_interval(1):h:z_interval(2);
    [y5] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y4(end,:)');
    
    % through quad 3
    z_interval = [params(10),params(10)+params(11)];
    h = h_quad;
    z6 = z_interval(1):h:z_interval(2);
    [y6] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y5(end,:)');
    
    % through drift 4
    z_interval = [params(10)+params(11),params(1)];
    h = h_drift;
    z7 = z_interval(1):h:z_interval(2);
    [y7] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y6(end,:)');
    
    if solenoid_flag
        % through solenoid
        z_interval = [params(1),params(1)+params(2)*0.50]; % 50% of the solenoid is fine
        h = h_sol;
        z8 = z_interval(1):h:z_interval(2);
        [y8] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(1), h, z_interval(2), y7(end,:)');
    else
        z8 = [];
        y8 = [];
    end
    
else
    
    h_drift = -h_drift;
    h_quad = -h_quad;
    h_sol = -h_sol;
    
    if solenoid_flag
        % through solenoid
        z_interval = [params(1),params(1)+params(2)*0.50]; % 50% of the solenoid is fine
        h = h_sol;
        z1 = z_interval(2):h:z_interval(1);
        [y1] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), init_cond);
    else
        z1 = [];
        y1 = init_cond';
    end
    
    % through drift 4
    z_interval = [params(10)+params(11),params(1)];
    h = h_drift;
    z2 = z_interval(2):h:z_interval(1);
    [y2] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y1(end,:)');
    
    % through quad 3
    z_interval = [params(10),params(10)+params(11)];
    h = h_quad;
    z3 = z_interval(2):h:z_interval(1);
    [y3] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y2(end,:)');
    
    % through drift 3
    z_interval = [params(7)+params(8),params(10)];
    h = h_drift;
    z4 = z_interval(2):h:z_interval(1);
    [y4] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y3(end,:)');
    
    % through quad 2
    z_interval = [params(7),params(7)+params(8)];
    h = h_quad;
    z5 = z_interval(2):h:z_interval(1);
    [y5] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y4(end,:)');
    
    % through drift 2
    z_interval = [params(4)+params(5),params(7)];
    h = h_drift;
    z6 = z_interval(2):h:z_interval(1);
    [y6] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y5(end,:)');
    
    % through quad 1
    z_interval = [params(4),params(4)+params(5)];
    h = h_quad;
    z7 = z_interval(2):h:z_interval(1);
    [y7] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y6(end,:)');
    
    % Drift distance up to first quad.
    z_interval = [0.0,params(4)];
    h = h_drift;
    z8 = z_interval(2):h:z_interval(1);
    [y8] = ode4(@(t,Y) odefcn(t,Y,params), z_interval(2), h, z_interval(1), y7(end,:)');
    
end

y = [y1;y2;y3;y4;y5;y6;y7;y8];
z = [z1,z2,z3,z4,z5,z6,z7,z8];

end

