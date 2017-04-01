//Function is composed by Su Dongcai on 2009/11/15
//If you have any suggestions, questions, and bug reports etc please feel free
//to contact me (suntree4152@gmail.com)
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "GraphSeg.h"

#ifndef mwSize
#define mwSize int
#endif

//const bool DEBUG=1;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
  mxArray *in_Img, *labeled_img, *threshold, *min_size, *radius, *knng, *knng_dist;
  double *p_inImg, *p_labeledImg, *pThreshold, *pMin_size, *pRadius, *pKnng, *pKnng_dist;
  int *pLabeledImg_out;
  const mwSize *dims;
  int width, height;
  int i;
  int numOvertex;
  bool bIsnnBased;
  if(nrhs==6)
  {
      bIsnnBased = 1;
  }
  else
  {
      bIsnnBased = 0;
  }
  //duplicate the input data
  in_Img = mxDuplicateArray(prhs[0]);
  threshold = mxDuplicateArray(prhs[1]);
  min_size = mxDuplicateArray(prhs[2]);
  radius = mxDuplicateArray(prhs[3]);
  if(bIsnnBased)
  {
      knng = mxDuplicateArray(prhs[4]);
      knng_dist = mxDuplicateArray(prhs[5]);
  }
  
  //access the input data by pointers:
  p_inImg = mxGetPr(in_Img);
  pThreshold = mxGetPr(threshold);
  pMin_size = mxGetPr(min_size);
  pRadius = mxGetPr(radius);
  if(bIsnnBased)
  {
      pKnng = mxGetPr(knng);
      pKnng_dist = mxGetPr(knng_dist);
  }
  
  
  //figure out the dimension of the input image
  dims = mxGetDimensions(prhs[0]);
  height = (int)dims[0];
  width = (int)dims[1];
  numOvertex = height*width;
  
  //construct the output data:
  labeled_img = plhs[0] = mxDuplicateArray(prhs[0]);
  p_labeledImg = mxGetPr(labeled_img);
  
  pLabeledImg_out = new int[numOvertex];
  if(bIsnnBased)
  {
      GraphSeg<double, int> segmentation(p_inImg, pLabeledImg_out, height, width, *pRadius, *pThreshold, *pMin_size, pKnng, pKnng_dist);
  }
  else
  {
      GraphSeg<double, int> segmentation(p_inImg, pLabeledImg_out, height, width, *pRadius, *pThreshold, *pMin_size);
  }
  
  
  //Give the output data:
  for(i=0; i<numOvertex; i++)
  {
      p_labeledImg[i] = (double)pLabeledImg_out[i];
  }
  
  //free up memory:
  
  delete [] pLabeledImg_out;
  
  mxDestroyArray(in_Img);
  mxDestroyArray(threshold);
  mxDestroyArray(min_size);
  mxDestroyArray(radius);
  if(bIsnnBased)
  {
      mxDestroyArray(knng);
      mxDestroyArray(knng_dist);
  }
  
  return ;
}