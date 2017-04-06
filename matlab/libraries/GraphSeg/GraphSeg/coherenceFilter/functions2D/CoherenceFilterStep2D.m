% This function CoherenceFilterStep2D is a MEX file which performs one 
% Anisotropic Diffusion  step on 2D gray/color image, using Weickerts 
% Rotation invariant diffusion scheme, and diffusion equation. 
%
% This mex file only uses CoherenceFilterStep2D.c and CoherenceFilterStep2D_functions.c
%
% Iout = CoherenceFilterStep2D(Iin, Options)
%
% inputs,
%   Iin : 2D gray/color image with datatype double
%   Options : Struct with filtering options (See CoherenceFilter.m)
%   
% outputs,
%   Iout : The anisotropic diffusion filtered image
%
% Written by D.Kroon University of Twente (September 2009)