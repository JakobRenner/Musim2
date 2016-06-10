function [ sv ] = von_Mises( sxx,syy,szz,sxy,syz,sxz )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% header for oofem .hom file:
% time exx eyy ezz eyz ezx exy sxx syy szz syz szx sxy      


 sv=sqrt(0.5*((sxx-syy).^2+(syy-szz).^2) + (sxx-szz).^2 + 6.*(syz.^2+sxz.^2+sxy.^2));

end

