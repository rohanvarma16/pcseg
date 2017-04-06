function [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=StructureTensor2DiffusionTensor3DWeickert(Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,Options)
% From Structure Tensor to Diffusion Tensor, a 3D implementation of the 2D 
% equations by Weickert
%
% [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=StructureTensor2DiffusionTensor3DWeickert(Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,Options)
% 
% Function is written by D.Kroon University of Twente (September 2009)

% Compute the eigenvectors and values of the structure tensors, v1, v2
% and v3, mu1, mu2 and mu3
[mu1,mu2,mu3,v3x,v3y,v3z,v2x,v2y,v2z,v1x,v1y,v1z]=EigenVectors3D(Jxx, Jxy, Jxz, Jyy, Jyz, Jzz);

[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=ConstructDiffusionTensor3D(v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,mu1,mu2,mu3,Options);

function  [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=ConstructDiffusionTensor3D(v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z,mu1,mu2,mu3,Options)
% Construct the edge preserving diffusion tensors D = [Dxx,Dxy,Dxz;Dxy,Dyy,Dyz;Dxz,Dyz,Dzz]

% Scaling of diffusion tensors
di=(mu1-mu3);
di((di<1e-15)&(di>-1e-15))=1e-15;
lambda1 = Options.alpha + (1 - Options.alpha)*exp(-Options.C./di.^(2*Options.m)); 
lambda2 = Options.alpha; 
lambda3 = Options.alpha;

% Construct the tensors
Dxx = lambda1.*v1x.^2   + lambda2.*v2x.^2   + lambda3.*v3x.^2;
Dyy = lambda1.*v1y.^2   + lambda2.*v2y.^2   + lambda3.*v3y.^2;
Dzz = lambda1.*v1z.^2   + lambda2.*v2z.^2   + lambda3.*v3z.^2;

Dxy = lambda1.*v1x.*v1y + lambda2.*v2x.*v2y + lambda3.*v3x.*v3y;
Dxz = lambda1.*v1x.*v1z + lambda2.*v2x.*v2z + lambda3.*v3x.*v3z;
Dyz = lambda1.*v1y.*v1z + lambda2.*v2y.*v2z + lambda3.*v3y.*v3z;

