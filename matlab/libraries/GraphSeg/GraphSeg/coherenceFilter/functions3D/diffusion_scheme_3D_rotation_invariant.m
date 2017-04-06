function u=diffusion_scheme_3D_rotation_invariant(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,dt)
% Most diffusion discretizations are not rotation-invariant. This is an 
% 3D standard discretization scheme which (implict) uses a 5x5x5 scheme, 
% based on the  2D scheme in "A Scheme for Coherence-Enhancing Diffusion 
% Filtering with Optimized Rotation invariance" by Joachim Weickert. 
% 
% Function is written by D.Kroon University of Twente (September 2009)

ux=derivatives(u,'x');
uy=derivatives(u,'y');
uz=derivatives(u,'z');

%% 3 : Calculate the flux components
j1 = Dxx .* ux + Dxy .*uy + Dxz .*uz;
j2 = Dxy .* ux + Dyy .*uy + Dyz .*uz;
j3 = Dxz .* ux + Dyz .*uy + Dzz .*uz;

%% 3.5 Set boundary flux to zero
j1(:,:,1)=0; j1(:,:,end)=0; j1(:,1,:)=0; j1(:,end,:)=0; j1(1,:,:)=0; j1(end,:,:)=0;
j2(:,:,1)=0; j2(:,:,end)=0; j2(:,1,:)=0; j2(:,end,:)=0; j2(1,:,:)=0; j2(end,:,:)=0;
j3(:,:,1)=0; j3(:,:,end)=0; j3(:,1,:)=0; j3(:,end,:)=0; j3(1,:,:)=0; j3(end,:,:)=0;


%% 4 : Calculate ... by means of the optimized derivative filters
du = derivatives(j1,'x')+derivatives(j2,'y')+derivatives(j3,'z');

%% 5 : Update in an explicit way.
u=u+du*dt;

	