function u=diffusion_scheme_3D_non_negativity(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,dt)
% This is a basic non_negativity discretization of the 3D diffusion
% filtering. 
%
% This implementation is based on fortran code from Achilles urangakis 
% and c-code from A Leith. Code and literature in software package 
% SPIDER (System for  Processing Image Data from Electron microscopy 
% and Related fields) 
% 
% Function is written by D.Kroon University of Twente (September 2009)

% Compute tensor-driven diffusion (as in [1] pp. 80-82)
[N1,N2,N3] = size(u);
px = [2:N1,N1]; nx = [1,1:N1-1];
py = [2:N2,N2]; ny = [1,1:N2-1];
pz = [2:N3,N3]; nz = [1,1:N3-1];

% Scales of image volume
HX=1; HY=1; HZ=1; 

% Some scaling constants
Rxx  = dt / (2.0 * HX * HX);
Ryy  = dt / (2.0 * HY * HY);
Rzz  = dt / (2.0 * HZ * HZ);
Rxy  = dt / (4.0 * HX * HY);
Rxz  = dt / (4.0 * HX * HZ);
Ryz  = dt / (4.0 * HY * HZ);
	  
% Stencil Weights
WE  = Rxx * (Dxx(px,:,:) + Dxx) - Rxy * (abs(Dxy(px,:,:)) + abs(Dxy)) - Rxz * (abs(Dxz(px,:,:)) + abs(Dxz));
WW  = Rxx * (Dxx(nx,:,:) + Dxx) - Rxy * (abs(Dxy(nx,:,:)) + abs(Dxy)) - Rxz * (abs(Dxz(nx,:,:)) + abs(Dxz));
WS  = Ryy * (Dyy(:,py,:) + Dyy) - Rxy * (abs(Dxy(:,py,:)) + abs(Dxy)) - Ryz * (abs(Dyz(:,py,:)) + abs(Dyz));
WN  = Ryy * (Dyy(:,ny,:) + Dyy) - Rxy * (abs(Dxy(:,ny,:)) + abs(Dxy)) - Ryz * (abs(Dyz(:,ny,:)) + abs(Dyz));
WB  = Rzz * (Dzz(:,:,nz) + Dzz) - Ryz * (abs(Dyz(:,:,nz)) + abs(Dyz)) - Rxz * (abs(Dxz(:,:,nz)) + abs(Dxz));
Wu  = Rzz * (Dzz(:,:,pz) + Dzz) - Ryz * (abs(Dyz(:,:,pz)) + abs(Dyz)) - Rxz * (abs(Dxz(:,:,pz)) + abs(Dxz));
WSE = Rxy * (   Dxy(px,py,:) + Dxy + abs(Dxy(px,py,:)) + abs(Dxy));
WNW = Rxy * (   Dxy(nx,ny,:) + Dxy + abs(Dxy(nx,ny,:)) + abs(Dxy));
WNE = Rxy * ( - Dxy(px,ny,:) - Dxy + abs(Dxy(px,ny,:)) + abs(Dxy));
WSW = Rxy * ( - Dxy(nx,py,:) - Dxy + abs(Dxy(nx,py,:)) + abs(Dxy));
WSu = Ryz * (   Dyz(:,py,pz) + Dyz + abs(Dyz(:,py,pz)) + abs(Dyz));
WNu = Ryz * ( - Dyz(:,ny,pz) - Dyz + abs(Dyz(:,ny,pz)) + abs(Dyz));
WEu = Rxz * (   Dxz(px,:,pz) + Dxz + abs(Dxz(px,:,pz)) + abs(Dxz));
WWu = Rxz * ( - Dxz(nx,:,pz) - Dxz + abs(Dxz(nx,:,pz)) + abs(Dxz));
WSB = Ryz * ( - Dyz(:,py,nz) - Dyz + abs(Dyz(:,py,nz)) + abs(Dyz));
WNB = Ryz * (   Dyz(:,ny,nz) + Dyz + abs(Dyz(:,ny,nz)) + abs(Dyz));
WEB = Rxz * ( - Dxz(px,:,nz) - Dxz + abs(Dxz(px,:,nz)) + abs(Dxz));
WWB = Rxz * (   Dxz(nx,:,nz) + Dxz + abs(Dxz(nx,:,nz)) + abs(Dxz));

% Perform the edge preserving diffusion filtering on the image volume
u = u ...
 + WE  .* (u(px,:,:)  - u)...
 + WW  .* (u(nx,:,:)  - u)...
 + WS  .* (u(:,py,:)  - u)...
 + WN  .* (u(:,ny,:)  - u)...
 + WB  .* (u(:,:,nz)  - u)...
 + Wu  .* (u(:,:,pz)  - u)...
 + WSE .* (u(px,py,:) - u)...
 + WNW .* (u(nx,ny,:) - u)...
 + WSW .* (u(nx,py,:) - u)...
 + WNE .* (u(px,ny,:) - u)...
 + WNB .* (u(:,ny,nz) - u)...
 + WNu .* (u(:,ny,pz) - u)...
 + WEB .* (u(px,:,nz) - u)...
 + WEu .* (u(px,:,pz) - u)...
 + WWB .* (u(nx,:,nz) - u)...
 + WWu .* (u(nx,:,pz) - u)...
 + WSB .* (u(:,py,nz) - u)...
 + WSu .* (u(:,py,pz) - u);
				 
