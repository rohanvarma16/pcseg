function [mu3,mu2,mu1,v3x,v3y,v3z,v2x,v2y,v2z,v1x,v1y,v1z]=EigenVectors3D(Dxx, Dxy, Dxz, Dyy, Dyz, Dzz)
% This function calculates the eigen values and vectors, of the 3D image hessian.
% 
% Function is written by D.Kroon University of Twente (September 2009)

v1x=zeros(size(Dxx),'single');
v1y=zeros(size(Dxx),'single');
v1z=zeros(size(Dxx),'single');
v2x=zeros(size(Dxx),'single');
v2y=zeros(size(Dxx),'single');
v2z=zeros(size(Dxx),'single');
v3x=zeros(size(Dxx),'single');
v3y=zeros(size(Dxx),'single');
v3z=zeros(size(Dxx),'single');
mu1=zeros(size(Dxx),'single');
mu2=zeros(size(Dxx),'single');
mu3=zeros(size(Dxx),'single');

for i=1:numel(Dxx)
    M=[Dxx(i) Dxy(i) Dxz(i); Dxy(i) Dyy(i) Dyz(i); Dxz(i) Dyz(i) Dzz(i)];
    [v,d]=eig(M);
    v=v';
    ev1=d(1,1); ev2=d(2,2); ev3=d(3,3);
    ev1a=abs(ev1); ev2a=abs(ev2); ev3a=abs(ev3);
    
    if((ev1a>=ev2a)&&(ev1a>ev3a))
        d=ev3; dt=ev3a;  dat=v(3,:);
        ev3=ev1; v(3,:)=v(1,:);
        ev1=d; ev1a=dt; v(1,:)=dat;
    elseif((ev2a>=ev1a)&&(ev2a>ev3a))
        d=ev3; dt=ev3a;  dat=v(3,:);
        ev3=ev1; v(3,:)=v(2,:);
        ev2=d; ev2a=dt; v(2,:)=dat;
    end
    if(ev1a>ev2a)
        d=ev2; dat=v(2,:);
        ev2=ev1; v(2,:)=v(1,:);
        ev1=d; v(1,:)=dat;
    end
    
    mu1(i)=ev1; mu2(i)=ev2; mu3(i)=ev3;
    v1x(i)=v(1,1); v1y(i)=v(1,2); v1z(i)=v(1,3);
    v2x(i)=v(2,1); v2y(i)=v(2,2); v2z(i)=v(2,3);
    v3x(i)=v(3,1); v3y(i)=v(3,2); v3z(i)=v(3,3);
end