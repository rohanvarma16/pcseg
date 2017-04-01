%effect: This script is for the test purpose of (Graph based segmentation algorithm)
%Composed by Su Dongcai on 2009/11/15
%If you have any suggestions, questions, and bug reports etc please feel free
%to contact me (suntree4152@gmail.com)

%Copyright (c) 2009, Su Dongcai
%All rights reserved.

%Acknowledgement:
%The Author want to thanks:
%[1]    Yan Zhang
%       2D/3D image segmentation toolbox
%       http://www.mathworks.com/matlabcentral/fileexchange/24998-2d3d-image-segmentation-toolbox
%For his "binaryHeap.h" class in cpp
%[2]    Dirk-Jan Kroon
%       Image Edge Enhancing Coherence Filter Toolbox
%       http://www.mathworks.com/matlabcentral/fileexchange/25449-image-edge-enhancing-coherence-filter-toolbox
%For the "Coherence Filter Toolbox" he provided
%And
%[3]    Luigi Giaccari
%       FAST K-NEAREST NEIGHBOURS SEARCH 3D VERSION
%       http://www.mathworks.com/matlabcentral/fileexchange/24607-fast-k-nearest-neighbours-search-3d-version
%For the fast implementation of K-NEAREST NEIGHBORS SEARCH in 3D.

% add all needed function paths
addpath .\coherenceFilter
addpath .\GLtree3DMex
%% Compile
fprintf('COMPILING:\n')
mex GraphSeg_mex.cpp
fprintf('\tGraphSeg_mex.cpp: mex succesfully completed.\n') 

mex .\GLtree3DMex\BuildGLTree.cpp
fprintf('\tBuildGLTree : mex succesfully completed.\n') 

mex .\GLtree3DMex\KNNSearch.cpp
fprintf('\tKNNSearch : mex succesfully completed.\n') 

mex .\GLtree3DMex\DeleteGLTree.cpp
fprintf('\tDeleteGLTree : mex succesfully completed.\n\n') 
%end of Complie#
%load an gray image:
load clown;
I_gray = X;
%smooth the image by coherence filter:
filted_I = CoherenceFilter(I_gray,struct('T',5,'rho',2,'Scheme','I', 'sigma', 1));
%adjacent neighborhood  model:
L = graphSeg(filted_I, 0.5, 50, 2, 0);
%k-nearest neighborhood model:
Lnn = graphSeg(filted_I, 0.5/sqrt(3), 50, 10, 1);
%display:
subplot(3, 1, 1), imshow(I_gray, []), title('original image');
subplot(3, 1, 2), imshow(label2rgb(L)), title('adjacent neighborhood based segmentation');
subplot(3, 1, 3), imshow(label2rgb(Lnn)), title('k nearest neighborhood based segmentation');