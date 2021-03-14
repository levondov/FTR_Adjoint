function [x,y] = solfield(pts,perc_trailing_edge)

x = linspace(0,1,pts);

x_left_t = 0:0.05:pi/2;
x_right_t = pi/2:0.05:pi;
y_left = cos(x_left_t+pi/2).^2;
y_right = cos(x_right_t+pi/2).^2;

x_left_edge_start = perc_trailing_edge;
x_right_edge_start = 1 - x_left_edge_start;

x_left = linspace(0,1,length(x_left_t))*x_left_edge_start;
x_right = linspace(0,1,length(x_left_t))*x_left_edge_start + x_right_edge_start;



% calculate y values
y = ones(pts,1);

for i = 1:pts
    if x(i) < x_left_edge_start % on the left trailing edge
        y(i) = interp1(x_left,y_left,x(i));
    end
    if x(i) > x_right_edge_start % on the right trailing edge
        y(i) = interp1(x_right,y_right,x(i));
    end
end




