function u=diffusion_scheme_2D_implicit(u,Dxx,Dxy,Dyy,dt)
% Diffusion scheme as introduced by Weickert  "Anisotropic Diffusion 
% in Image Processing", Thesis 1996.
%
% This Implementation is based on Pascal Getreuer, "Weickert's
% coherence-enhancing filter" which is also on Mathworks.
% 
% Function is written by D.Kroon University of Twente (September 2009)

% Compute tensor-driven diffusion (as in [1] pp. 80-82)
[N1,N2,N3] = size(u);
id = [2:N1,N1]; iu = [1,1:N1-1];
ir = [2:N2,N2]; il = [1,1:N2-1];

% Diagonal neighbor intensity updates (amount of greyvalue flow from
% neighbor) 

% Page 95
% This formula is equal to ( |b       | - b           +   |b    | - b    ) / ( 4*h1*h2 )
%                            | i-1,j+1|    i-1,j+1        | i,j |    i,j
Dul = max(0, Dxy(iu,il)) + max(0,Dxy); 
Ddr = max(0, Dxy(id,ir)) + max(0,Dxy);
Ddl = max(0,-Dxy(id,il)) + max(0,-Dxy);
Dur = max(0,-Dxy(iu,ir)) + max(0,-Dxy);

% Normal neighbor intensity updates 
Dml = (Dyy(:,il) + Dyy) - (abs(Dxy(:,il)) + abs(Dxy));
Dmr = (Dyy(:,ir) + Dyy) - (abs(Dxy(:,ir)) + abs(Dxy));
Duc = (Dxx(iu,:) + Dxx) - (abs(Dxy(iu,:)) + abs(Dxy));
Ddc = (Dxx(id,:) + Dxx) - (abs(Dxy(id,:)) + abs(Dxy));

% Normalization term to preserve region average grey value
Dmc = -(Dyy(:,il) + 2*Dyy + Dyy(:,ir)) + ... 
      -(Dxx(iu,:) + 2*Dxx + Dxx(id,:)) + ...  
      -(max(0, Dxy(iu,il)) + max(0,-Dxy(iu,ir))) + ...
      -(max(0,- Dxy(id,il)) + max(0,Dxy(id,ir))) + ...
      (abs(Dxy(:,il)) + abs(Dxy(:,ir)) + abs(Dxy(id,:)) + abs(Dxy(iu,:)) + 2*abs(Dxy));

% Calculate Diffusion update 
du=zeros(size(u));
for Channel = 1:N3
    du(:,:,Channel)=0.5*dt*(Dul.*u(iu,il,Channel) + ...
                            Dml.*u(:,il,Channel)  + ...
                            Ddl.*u(id,il,Channel) + ...
                            Duc.*u(iu,:,Channel)  + ...            
                            Ddc.*u(id,:,Channel)  + ...
                            Dur.*u(iu,ir,Channel) + ...
                            Dmr.*u(:,ir,Channel)  + ...
                            Ddr.*u(id,ir,Channel));
end

% Perform the edge preserving diffusion filtering on the image
for Channel = 1:N3
    u(:,:,Channel) = (u(:,:,Channel) + du(:,:,Channel)) ./ (1-dt*0.5*Dmc);
end

