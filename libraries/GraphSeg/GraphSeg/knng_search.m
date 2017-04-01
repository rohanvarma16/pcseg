function [knng, knng_dist] = knng_search(img, numOnns)
%Inputs:
%       img: the normalized double img in [0, 1]
%       numOnns: the number of the nearest neighbor points of a pixel.
%Ouputs:
%       knng: the idx matrix of the knn of each pixel
%       knng_dist: the corresponding distance of knng

%Composed by Su Dongcai on 2009/11/15
%If you have any suggestions, questions, and bug reports etc please feel free
%to contact me (suntree4152@gmail.com)

%Copyright (c) 2009, Su Dongcai
%All rights reserved.

height = size(img, 1);
width = size(img, 2);
numOpixels = height*width;
Ref_set = zeros(numOpixels, 3);
Ref_idx=1;
for x=1:width
    for y=1:height
        gray = img(y, x);
        Ref_set(Ref_idx, :)=[x, y, gray];
        Ref_idx = Ref_idx+1;
    end
end

Ref_set(:, 1)=Ref_set(:, 1)/width;
Ref_set(:, 2)=Ref_set(:, 2)/height;
Query_set = Ref_set;
%begin knn search:
pTree = BuildGLTree(Ref_set);
[knng, knng_dist] = KNNSearch(Ref_set, Query_set, pTree, numOnns);
%#debug:
%min(knng_dist)
knng = knng-1;
DeleteGLTree(pTree);