% This function CoherenceFilterStep3D is a MEX file which performs one 
% Anisotropic Diffusion  step on 3D image volume, and is an extention of
% Weickerts 2D Rotation invariant diffusion scheme, and 2D diffusion equation. 
%
% This mex file only uses CoherenceFilterStep3D.c and CoherenceFilterStep3D_functions.c
%
% Iout = CoherenceFilterStep2D(Iin, Options)
%
% inputs,
%   Iin : 3D image with datatype single
%   Options : Struct with filtering options (See CoherenceFilter.m)
%   
% outputs,
%   Iout : The anisotropic diffusion filtered image
%
% Written by D.Kroon University of Twente (September 2009)