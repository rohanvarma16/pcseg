function L = graphSeg(img, threshold, min_size, nRadius, model)
%Input:
%       img: the gray image
%       threshold: larger prefer larger segmented area
%       min_size: the minimum size of segmentation component
%       nRadius: the radius of neighbourhood of a pixel if in model (0)
%       nRadius: the number of nearest neighborhood of a pixel if in model
%       (1)
%       model: 0-->adjacent neighborhood based
%              1-->k nearest neighborhood based
%       Note: the precisely meaning of above parameters please refer to [1]
%Output:
%       L: the labeled image, differente area is labeled by different
%       number
%Examples:
% %load an gray image:
% load clown;
% I_gray = X;
% %smooth the image by coherence filter:
% filted_I = CoherenceFilter(I_gray,struct('T',5,'rho',2,'Scheme','I', 'sigma', 1));
% %adjacent neighborhood  model:
% L = graphSeg(filted_I, 0.5, 50, 2, 0);
% %k-nearest neighborhood model:
% Lnn = graphSeg(filted_I, 0.5/sqrt(3), 50, 10, 1);
% %display:
% subplot(3, 1, 1), imshow(I_gray, []), title('original image');
% subplot(3, 1, 2), imshow(label2rgb(L)), title('adjacent neighborhood based segmentation');
% subplot(3, 1, 3), imshow(label2rgb(Lnn)), title('k nearest neighborhood based segmentation');

%Reference:
%[1] Pedro F. Felzenszwalb and Daniel P. Huttenlocher 
%    International Journal of Computer Vision, Volume 59, Number 2, September 2004
%[2] http://people.cs.uchicago.edu/~pff/segment/

%Composed by Su Dongcai on 2009/11/15
%If you have any suggestions, questions, and bug reports etc please feel free
%to contact me (suntree4152@gmail.com)

%Copyright (c) 2009, Su Dongcai
%All rights reserved.
img=im2double(img);
img = img/max(img(:));
if model==0
    L = GraphSeg_mex(img, threshold, min_size, nRadius);
elseif model==1
    [knng, knng_dist] = knng_search(img, nRadius);
    L = GraphSeg_mex(img, threshold, min_size, nRadius, knng, knng_dist);
else
    L = [];
end