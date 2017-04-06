% BuildGLTree construct a GL-Tree from a 2D point cloud
%
% SYNTAX
% ptrtree=BuildGLTree(px,py);
%
% INPUT PARAMETERS
%   px,py: [Nx1] double vectors of the x and y coordinates of points.
%     
%
% OUTPUT PARAMETERS
%   ptrtree: a pointer to the created data structure
%
%
% GENERAL INFORMATIONS
%
%     - GLTree is an exact method no approximation is done. If you find a
%      different value from the expected this means you found a bug so please
%      send a report to the author.
%     - GLTree works on double precision so these must be double vectors.
%     - GLTree do not support 3D points.
%     - The Data structure will be computed in linear time with the number of points.
%     - GLTree is faster on uniformly random data. On sparse ones should work
%      properly but may be  slower
%
%
%For question, suggestion, bug reports
%giaccariluigi@msn.com
%
% Visit my website:
% http://giaccariluigi.altervista.org/blog/
%
%Author : Luigi Giaccari
%Last Update: 7/12/2008
%Created : 10/10/2008