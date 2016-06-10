function [ f_n ] = normal_force( Y, nu, A, R, r, v )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Y = Y * 10e5
R_eff = R(1) * R(2) / (R(1) + R(2))
u = norm( r(:,1)-r(:,2) )

rho = (2 * Y * R_eff^0.5) / (3 *(1-nu^2))

xi = R(1) + R(2) - u

dt_xi = - ( (r(1,1)-r(2,1))*(v(1,1)-v(2,1))/u + (r(1,2)-r(2,2))*(v(1,2)-v(2,2))/u ) 

f_n = rho * xi^(1.5) + 1.5 * rho * A * dt_xi * xi^0.5

end

function [ pos ] = contact_search(filename, radius)
fileID = fopen(filename,'r');
formatSpec = '%f';
pos = fscanf(fileID,formatSpec)
