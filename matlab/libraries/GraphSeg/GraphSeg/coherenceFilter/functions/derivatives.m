function D=derivatives(I,option)
% Sobel like derivatives with Scharr rotation invariance stencil notations.
%
% D=derivatives(I,option)
%
% inputs,
%   I : Input image 2D or 3D
%   option : Derivative direction 'x', 'y' or 'z'
%
% outputs,
%   D : Derivative of input image
%
% Written by D.Kroon University of Twente (September 2009)

if((ndims(I)<2)||ndims(I)>3)
    error('derivatives:unknowdim','Only 2D and 3D supported');
end
is3d=(size(I,3)>3);

m=1; sigma=1;
switch lower(option)
    case 'x'
        if(is3d)
            if(m==1)
                Stencil=zeros(3,3,3);
                Stencil(:,:,1)=[9 30 9; 0 0 0; -9 -30 -9];
                Stencil(:,:,2)=[30 100 30; 0 0 0; -30 -100 -30];
                Stencil(:,:,3)=[9 30 9; 0 0 0; -9 -30 -9];
                Stencil=Stencil/512;
            elseif(m==2)
                b=floor(-3*sigma):ceil(3*sigma); [x,y,z]=ndgrid(b,b,b);
                Stencil=-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));
            elseif(m==3)
                Stencil=[1 0 -1];
            end
        else
            if(m==1)
                Stencil= [3 10 3; 0 0 0; -3 -10 -3]/32;
            elseif(m==2)
                b=floor(-3*sigma):ceil(3*sigma); [x,y]=ndgrid(b,b);
                Stencil=-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
            elseif(m==3)
                Stencil=[1 0 -1];
            end
        end
    case 'y'
        if(is3d)
            if(m==1)
                Stencil=zeros(3,3,3);
                Stencil(:,:,1)=[9 0 -9; 30 0 -30; 9 0 -9];
                Stencil(:,:,2)=[30  0 -30; 100 0 -100; 30  0 -30];
                Stencil(:,:,3)=[9 0 -9; 30 0 -30; 9 0 -9];
                Stencil=Stencil/512;
            elseif(m==2)
                b=floor(-3*sigma):ceil(3*sigma); [x,y,z]=ndgrid(b,b,b);
                Stencil=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));
            elseif(m==3)
                Stencil=[1;0;-1];
            end
        else
            if(m==1)
                Stencil= [3  0 -3; 10 0 -10; 3  0 -3]/32;
            elseif(m==2)
                b=floor(-3*sigma):ceil(3*sigma); [x,y]=ndgrid(b,b);
                Stencil=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
            elseif(m==3)
                Stencil=[1;0;-1];
            end
        end
    case 'z'
        if(is3d)
            if(m==1)
                Stencil=zeros(3,3,3);
                Stencil(:,:,1)=[9  30  9; 30 100 30; 9  30  9];
                Stencil(:,:,2)=[0 0 0; 0 0 0; 0 0 0];
                Stencil(:,:,3)=[-9  -30  -9; -30 -100 -30; -9  -30  -9];
                Stencil=Stencil/512;
            elseif(m==2)
                b=floor(-3*sigma):ceil(3*sigma); [x,y,z]=ndgrid(b,b,b);
                Stencil=-(z./(2*pi*sigma^4)).*exp(-(x.^2+y.^2+z.^2)/(2*sigma^2));
            elseif(m==3)
                Stencil=zeros(1,1,3);
                Stencil(:,:,1)=1; Stencil(:,:,3)=-1;
            end
            
        else
            error('derivatives:unknowopt','2D thus no z derivative');
        end
    otherwise
        error('derivatives:unknowopt','Unknow option found');
end
D=imfilter(I,Stencil,'conv','same','replicate');

