function [M] = calcmoment(fxy,p,n)
% image moment calculation
%
% INPUTS:
%   fxy - a Nx2 array of particles column 1 is x pos. column 2 is y pos.
%   p - x coefficient
%   n - y coefficient
%
% OUTPUTS:
%   M - moment
%
%
% M_pn = \sum_x \sum_y x^p y^n I(x,y)
%
N = length(fxy(:,1));

if p == 1 && n == 0
    M = mean(fxy(:,1));
elseif p == 0 && n == 1
    M = mean(fxy(:,2));
elseif p == 2 && n == 0
    M = 0;
    xbar=mean(fxy(:,1));
    for i = 1:N
       M = M + abs(fxy(i,1)-xbar)^2;
    end
    M = M/(N);
elseif p == 0 && n == 2
    M = 0;
    ybar=mean(fxy(:,2));
    for i = 1:N
       M = M + abs(fxy(i,2)-ybar)^2;
    end
    M = M/(N);
elseif p == 1 && n == 1
    M = 0;
    xbar=mean(fxy(:,1));
    ybar=mean(fxy(:,2));
    for i = 1:N
        M = M + (fxy(i,1)-xbar)*(fxy(i,2)-ybar);
    end
    M = M/N;
end

end


