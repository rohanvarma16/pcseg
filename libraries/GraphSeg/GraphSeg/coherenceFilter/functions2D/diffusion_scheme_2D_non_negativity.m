function u=diffusion_scheme_2D_non_negativity(u,Dxx,Dxy,Dyy,dt)
% The Basic non_negativity diffusion equation. (Can be found in "A Scheme for
% Coherence-Enhancing Diffusion Filtering with Optimized Rotation
% Invariance" by Joachim Weickert
% 
% Function is written by D.Kroon University of Twente (September 2009)

% Make positive and negative indices
[N1,N2,N3] = size(u);
px = [2:N1,N1]; nx = [1,1:N1-1];
py = [2:N2,N2]; ny = [1,1:N2-1];

% In literature a,b and c are used as variables 
a=Dxx;b=Dxy;c=Dyy;

% Stencil Weights
A1 = (1/4)*((abs(b(nx, py))-b(nx,py)) + (abs(b)-b));
A2 = (1/2)*( (c(:, py)+c) -(abs(b(:,py))+abs(b)));
A3 = (1/4)*((abs(b(px, py))+b(px,py)) + (abs(b)+b));
A4 = (1/2)*( (a(nx,: )+a) -(abs(b(nx,:))+abs(b)));
%A5 =  - (a(nx,:) + 2*a+a(px,:)) ...
%	  - (abs(b(nx,py))-b(nx,py)+abs(b(px,py))-b(px,py)) ...
%	  - (abs(b(nx,ny))-b(nx,ny)+abs(b(px,ny))-b(px,ny)) ...
%	  + (abs(b(nx,:))+abs(b(px,:))+abs(b(:,ny))+abs(b(:,py))+2*abs(b)) ...
%	  - (c(:,ny) + 2*c+c(:,py));
A6 = (1/2)*( (a(px,: )+a) -(abs(b(px,:))+abs(b)));
A7 = (1/4)*((abs(b(nx, ny))+b(nx,ny)) + (abs(b)+b));
A8 = (1/2)*( (c(:, ny)+c) -(abs(b(:,ny))+abs(b)));
A9 = (1/4)*((abs(b(px, ny))-b(px,ny)) + (abs(b)-b));

% Do the diffusion
for Channel = 1:N3
  u(:,:,Channel)=  u(:,:,Channel)+dt*(A1.*(u(nx,py,Channel) -u(:,:,Channel))+ ...
                                   A2.*(u(:, py,Channel) -u(:,:,Channel))+ ...
                                   A3.*(u(px,py,Channel) -u(:,:,Channel))+ ...
                                   A4.*(u(nx,:,Channel)  -u(:,:,Channel))+ ...      
                                   A6.*(u(px,:,Channel)  -u(:,:,Channel))+ ...
                                   A7.*(u(nx,ny,Channel) -u(:,:,Channel))+ ...
                                   A8.*(u(:, ny,Channel) -u(:,:,Channel))+ ...
                                   A9.*(u(px,ny,Channel) -u(:,:,Channel)));
end
