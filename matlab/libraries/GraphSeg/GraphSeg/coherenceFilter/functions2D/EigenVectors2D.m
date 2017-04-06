function [mu1,mu2,v1x,v1y,v2x,v2y]=EigenVectors2D(Jxx,Jxy,Jyy)
% This function computes the eigenvectors and eigen values of the 2D image
% Hessian
%
% [mu1,mu2,v1x,v1y,v2x,v2y]=EigenVectors2D(Jxx,Jxy,Jyy)
%
% inputs, 
%   Jxx, Jxy and Jyy : Matrices with the values of the Hessian tensors
% 
% outputs,
%   mu1, mu2 : Matrices with eigen values
%   v1x, v1y, v2x, v2y : Matrices with the eigen vectors
% 
% Function is written by D.Kroon University of Twente (September 2009)

% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Jxx - Jyy).^2 + 4*Jxy.^2);
v2x = 2*Jxy; v2y = Jyy - Jxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y; v1y = v2x;

% Compute the eigenvalues
mu1 = 0.5*(Jxx + Jyy + tmp);
mu2 = 0.5*(Jxx + Jyy - tmp);

