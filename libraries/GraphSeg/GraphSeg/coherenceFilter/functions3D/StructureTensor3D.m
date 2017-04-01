function [Jxx, Jxy, Jxz, Jyy, Jyz, Jzz]=StructureTensor3D(ux,uy,uz, rho)
% This function calculates the 3D Structure Tensor
% Jp ( grad(u) ) = conv ( Kp  , grad(u) * grad(u)^T )
% 
% Function is written by D.Kroon University of Twente (September 2009)

% J(grad u_sigma)
Jxx = ux.^2;
Jxy = ux.*uy;
Jxz = ux.*uz;
Jyy = uy.^2;
Jyz = uy.*uz;
Jzz = uz.^2;

% Do the gaussian smoothing
Jxx = imgaussian(Jxx,rho,4*rho);
Jxy = imgaussian(Jxy,rho,4*rho);
Jxz = imgaussian(Jxz,rho,4*rho);
Jyy = imgaussian(Jyy,rho,4*rho);
Jyz = imgaussian(Jyz,rho,4*rho);
Jzz = imgaussian(Jzz,rho,4*rho);
