function u = CoherenceFilter(u,Options)
% This function COHERENCEFILTER will perform Anisotropic Diffusion of a
% 2D gray/color image or 3D image volume, Which will reduce the noise in
% an image while preserving the region edges, and will smooth along
% the image edges removing gaps due to noise.
%
% Don't forget to compile the c-code by executing compile_c_files.m
% 
% Iout = CoherenceFilter(Iin, Options)
%
% inputs,
%   Iin : 2D gray/color image or 3D image volume. Use double datatype in 2D
%              and single data type in 3D. Range of image data must be 
%              approximately [0 1], for default constants.
%   Options : Struct with filtering options
%   
% outputs,
%   Iout : The anisotropic diffusion filtered image
%
% Options,
%   Options.Scheme :  The numerical diffusion scheme used
%                     'R', Rotation Invariant, Standard Discretization 
%                          (implicit) 5x5 kernel (Default)
%                     'I', Implicit Discretization (only works in 2D)
%                     'S', Standard Discretization
%                     'N', Non-negativity Discretization
%   Options.T  :      The total diffusion time (default 5)
%   Options.dt :      Diffusion time stepsize, in case of scheme R or I
%                     defaults to 1, in case of scheme S or N defaults to
%                     0.1. 
%   Options.sigma :   Sigma of gaussian smoothing before calculation of the
%                     image Hessian, default 1.                   
%   Options.rho :     Rho gives the sigma of the Gaussian smoothing of the 
%                     Hessian, default 1.
%   Options.verbose : Show information about the filtering, values :
%                     'none', 'iter' (default) , 'full'
%
% Constants which determine the amplitude of the diffusion smoothing in 
% Weickert equation
%   Options.C :     Default 1e-10
%   Options.m :     Default 1
%   Options.alpha : Default 0.001
%                   
% The basis of the method used is the one introduced by Weickert:
%   1, Calculate Hessian from every pixel of the gaussian smoothed input image
%   2, Gaussian Smooth the Hessian, and calculate its eigenvectors and values
%      (Image edges give large eigenvalues, and the eigenvectors corresponding
%         to those large eigenvalues describe the direction of the edge)
%   3, The eigenvectors are used as diffusion tensor directions. The 
%      amplitude of the diffusion in those 3 directions is determined
%      by equations below.
%   4, An Finite Difference scheme is used to do the diffusion
%   5, Back to step 1, till a certain diffusion time is reached.
%
% Weickert equation 2D:
%    lambda1 = alpha + (1 - alpha)*exp(-C/(mu1-mu2).^(2*m)); 
%    lambda2 = alpha;
% Weickert extended to 3D:
%    lambda1 = alpha + (1 - alpha)*exp(-C/(mu1-mu3).^(2*m)); 
%    lambda2 = alpha;
%    lambda3 = alpha;
%  (with mu1 the largest eigenvalue and mu3 the smallest)
%
% Notes:
% - If the time step is choosing to large the scheme becomes unstable, this
%   can be seen by setting verbose to 'full'. The image variance has to
%   decrease every itteration if the scheme is stable.
% - Weickert's equation can be found in several code files. But this is
%   only one of the possible diffusion equations, you can for instance
%   do plane smoothing instead of line smoothing by edditing "lambda2 = alpha"
%   to "lambda2 = alpha + (1 - alpha)*exp(-C/(mu2-mu3).^(2*m)); ", or
%   use the equation of Siham Tabik.
%
% Literature used:
%  - Weickert : "A Scheme for Coherence-Enhancing Diffusion Filtering
%                   with Optimized Rotation Invariance"
%  - Weickert : "Anisotropic Diffusion in Image Processing", Thesis 1996
%  - Laura Fritz : "Diffusion-Based Applications for Interactive Medical
%                   Image Segmentation"
%  - Siham Tabik, et al. : "Multiprocessing of Anisotropic Nonlinear
%                          Diffusion for filtering 3D image"
%  
% example 2d,
%   I = im2double(imread('images/sync_noise.png'));
%   JR = CoherenceFilter(I,struct('T',4,'rho',10,'Scheme','R'));
%   JI = CoherenceFilter(I,struct('T',4,'rho',10,'Scheme','I'));
%   JS = CoherenceFilter(I,struct('T',4,'rho',10,'Scheme','N'));
%   figure, 
%   subplot(2,2,1), imshow(I), title('Before Filtering');
%   subplot(2,2,2), imshow(JR), title('Rotation Invariant Scheme');
%   subplot(2,2,3), imshow(JI), title('Implicit Scheme');
%   subplot(2,2,4), imshow(JI), title('Non Negative Scheme');
%
% example 3d,
%	% First compile the c-code by executing compile_c_files.m
%   load('images/sphere');
%   showcs3(V);
%   JR = CoherenceFilter(V,struct('T',50,'dt',2,'Scheme','R'));
%   showcs3(JR);
%
% Written by D.Kroon University of Twente (September 2009)

% add all needed function paths
try
    functionname='CoherenceFilter.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions2D'])
    addpath([functiondir '/functions3D'])
    addpath([functiondir '/functions'])
catch me
    disp(me.message);
end

% Default parameters
defaultoptions=struct('T',2,'dt',[],'sigma', 1, 'rho', 1, 'TensorType', 1, 'C', 1e-10, 'm',1,'alpha',0.001,'C2',0.3,'m2',8,'alpha2',1,'RealDerivatives',false,'Scheme','R','verbose','iter');

if(~exist('Options','var')),
    Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))),
        warning('CoherenceFilter:unknownoption','unknown options found');
    end
end

if(isempty(Options.dt))
    switch lower(Options.Scheme)
      case 'r', Options.dt=1;
      case 'i', Options.dt=1;
      case 's', Options.dt=0.1;
      case 'n', Options.dt=0.1;
      otherwise
        error('CoherenceFilter:unknownoption','unknown scheme');
    end
end
    

% Initialization
dt_max = Options.dt; t = 0;

% In case of 3D use single precision to save memory
if(size(u,3)<4), u=double(u); else u=single(u); end 

% Process time
%process_time=tic;
tic;
% Show information 
switch lower(Options.verbose(1))
case 'i'
    disp('Diffusion time   Sec. Elapsed');
case 'f'
    disp('Diffusion time   Sec. Elapsed   Image mean    Image variance');
end        


% Anisotropic diffusion main loop
while (t < (Options.T-0.001))
    % Update time, adjust last time step to exactly finish at the wanted
    % diffusion time
    Options.dt = min(dt_max,Options.T-t); t = t + Options.dt;
    tn=toc;
    switch lower(Options.verbose(1))
    case 'n'
    case 'i'
        s=sprintf('    %5.0f        %5.0f    ',t,round(tn)); disp(s);
    case 'f'
        s=sprintf('    %5.0f        %5.0f      %13.6g    %13.6g ',t,round(tn), mean(u(:)), var(u(:))); disp(s);
    
    end        
   
    if(size(u,3)<4) % Check if 2D or 3D
        % Do a diffusion step
        if(strcmpi(Options.Scheme,'R@'))
            u=CoherenceFilterStep2D(u,Options);
        else
            u=Anisotropic_step2D(u,Options);
        end
    else
        % Do a diffusion step
        if(strcmpi(Options.Scheme,'R'))
            u=CoherenceFilterStep3D(u,Options);
        else
            u=Anisotropic_step3D(u,Options);
        end
    end
end

function u=Anisotropic_step2D(u,Options)
% Perform tensor-driven diffusion filtering update

% Gaussian smooth the image, for better gradients
usigma=imgaussian(u,Options.sigma,4*Options.sigma);


% Calculate the gradients
switch lower(Options.Scheme)
  case 'r'
    ux=derivatives(usigma,'x'); uy=derivatives(usigma,'y');
  case 'i'
    ux=derivatives(usigma,'x'); uy=derivatives(usigma,'y');
  case 's'
    [uy,ux]=gradient(usigma);
  case 'n'
    [uy,ux]=gradient(usigma);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end


% Compute the 2D structure tensors J of the image
[Jxx, Jxy, Jyy] = StructureTensor2D(ux,uy,Options.rho);

% Compute the eigenvectors and values of the strucure tensors, v1 and v2, mu1 and mu2
[mu1,mu2,v1x,v1y,v2x,v2y]=EigenVectors2D(Jxx,Jxy,Jyy);

% Construct the edge preserving diffusion tensors D = [Dxx,Dxy;Dxy,Dyy]
[Dxx,Dxy,Dyy]=ConstructDiffusionTensor2D(mu1,mu2,v1x,v1y,v2x,v2y,Options);

% Do the image diffusion
switch lower(Options.Scheme)
  case 'r'
      u=diffusion_scheme_2D_rotation_invariant(u,Dxx,Dxy,Dyy,Options.dt);
  case 'i'
      u=diffusion_scheme_2D_implicit(u,Dxx,Dxy,Dyy,Options.dt);
  case 's'
      u=diffusion_scheme_2D_standard(u,Dxx,Dxy,Dyy,Options.dt);
  case 'n'
      u=diffusion_scheme_2D_non_negativity(u,Dxx,Dxy,Dyy,Options.dt);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end

function u=Anisotropic_step3D(u,Options)
% Perform tensor-driven diffusion filtering update

% Gaussian smooth the image, for better gradients
usigma=imgaussian(u,Options.sigma,4*Options.sigma);

% Calculate the gradients
ux=derivatives(usigma,'x');
uy=derivatives(usigma,'y');
uz=derivatives(usigma,'z');

% Compute the 3D structure tensors J of the image
[Jxx, Jxy, Jxz, Jyy, Jyz, Jzz] = StructureTensor3D(ux,uy,uz, Options.rho);

% Free memory
clear ux; clear uy; clear uz;

% Compute the eigenvectors and eigenvalues of the hessian and directly
% use the equation of Weickert to convert them to diffusion tensors
[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=StructureTensor2DiffusionTensor3DWeickert(Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,Options); 

% Free memory
clear J*;

% Do the image diffusion
switch lower(Options.Scheme)
  case 'r'
      u=diffusion_scheme_3D_rotation_invariant(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 'i'
      u=diffusion_scheme_3D_implicit(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 's'
      u=diffusion_scheme_3D_standard(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 'n'
      u=diffusion_scheme_3D_non_negativity(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end




        