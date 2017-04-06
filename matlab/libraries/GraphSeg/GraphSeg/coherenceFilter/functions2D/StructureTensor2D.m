function [Jxx, Jxy, Jyy]=StructureTensor2D(ux,uy,rho)
% This function calculates the 2D
% Jp ( grad(u) ) = conv ( Kp  , grad(u) * grad(u)^T ) 
%
% The 2x2 JP matrix is called structure tensor (second-moment matrix,
% scatter matrix, Forstner interest operator)
% J. Weickert "Coherence-Enhancing Shock filters"
% 
% Function is written by D.Kroon University of Twente (September 2009)

% J(grad u_sigma)
Jxx = sum(ux.^2,3)/size(ux,3);
Jxy = sum(ux.*uy,3)/size(ux,3);
Jyy = sum(uy.^2,3)/size(ux,3);

% Do the gaussian smoothing of the structure tensor
Jxx = imgaussian(Jxx,rho,6*rho);
Jxy = imgaussian(Jxy,rho,6*rho);
Jyy = imgaussian(Jyy,rho,6*rho);