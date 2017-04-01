function [Dxx,Dxy,Dyy]=ConstructDiffusionTensor2D(mu1,mu2,v1x,v1y,v2x,v2y,Options)
% Construct the edge preserving diffusion tensor D = [a,b;b,c] such as
% introduced by Weickert.
%
% http://www.mia.uni-saarland.de/weickert/Papers/book.pdf, pp 127-128
% 
% Function is written by D.Kroon University of Twente (September 2009)

di=(mu1-mu2);
di((di<1e-15)&(di>-1e-15))=1e-15;
% Implicit if mu1 == mu2 then lambda1=alpha
lambda1 = Options.alpha + (1 - Options.alpha)*exp(-Options.C./(di).^(2*Options.m)); 
lambda2 = Options.alpha;
 
Dxx = lambda1.*v1x.^2   + lambda2.*v2x.^2;
Dxy = lambda1.*v1x.*v1y + lambda2.*v2x.*v2y;
Dyy = lambda1.*v1y.^2   + lambda2.*v2y.^2;




