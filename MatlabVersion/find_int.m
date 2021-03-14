function [int_val] = find_int(z,y_adj,y,O_opt,N_opt)
%FIND_INT Summary of this function goes here

global norm_flag
% calculate adjoint integral broken up into 4 pieces
int1 = zeros(1,length(z));
int2 = zeros(1,length(z));
int3 = zeros(1,length(z));
int4 = zeros(1,length(z));
for i = 1:length(z)
   int1(i) = y_adj(i,4:6)*(O_opt{i}*y(i,1:3)');
   int2(i) = y_adj(i,10)*(y(i,1:3)*N_opt{i});
   int3(i) = -1*y_adj(i,1:3)*(O_opt{i}*y(i,4:6)');
   int4(i) = -1*y_adj(i,1:3)*N_opt{i}*y(i,10);
end

global k0

int_val = trapz(z,int1+int2+int3+int4);

if norm_flag
    int_val = int_val / ((k0^2) * (y(1,1)^2));
end

end

