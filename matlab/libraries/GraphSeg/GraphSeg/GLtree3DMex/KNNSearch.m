% KNNSearch query a GL-tree for k nearest neighbor(kNN)
%
% SYNTAX
%
% [kNNG]=KNNSearch(,qp,ptrtree,k);       short
% [kNNG,Dist]=KNNSearch(p,qp,ptrtree,k);  long
%
% INPUT PARAMETERS
% 
%       p: [Nx3] double array coordinates of reference points
% 
%       qp: [Nqx3] double array coordinates of query points
%
%       ptrtree: a pointer to the previously constructed  GLtree.Warning
%                if the pointer is uncorrect it will cause a crash, there is
%                no way to check this in the mex routine, you have to check
%                yourself in your script.
%
%       k: number of neighbors
%
% OUTPUT PARAMETERS
%
%      kNNG: [Nqxk] array, each rows contains the kNN indexes
%            So in row one there are kNN to first query
%           point, in row two to the second etc...
% 
%      Dist: [Nqxk] array, Facultative output, each rows contains the
%                   distance values of the  found kNN.
%         
%
%
%
% GENERAL INFORMATIONS
%
%         -This function is faster if all query points are given once
%         instead of looping and pass one point each loop.
%
%
%  For question, suggestion, bug reports
%  giaccariluigi@msn.com
% 
% Visit my website:
% http://giaccariluigi.altervista.org/blog/
%
%  Author : Luigi Giaccari
%  Last Update: 2/1/2009
%  Created : 8/8/2008
