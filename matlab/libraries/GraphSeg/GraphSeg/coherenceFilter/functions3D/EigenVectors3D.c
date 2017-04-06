#include "mex.h"
#include "math.h"
#include "EigenDecomposition3.c"

#ifndef mwSize
#define mwSize int
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    float *Dxx, *Dxy, *Dxz, *Dyy, *Dyz, *Dzz;
    float *Deiga, *Deigb, *Deigc;
	float *Dvecxa, *Dvecya, *Dvecza;
	float *Dvecxb, *Dvecyb, *Dveczb;
	float *Dvecxc, *Dvecyc, *Dveczc;

    mwSize output_dims[2]={1, 3};
    double Ma[3][3];
    double Davec[3][3];
    double Daeig[3];
    
    /* Loop variable */
    int i;
    
    /* Size of input */
    const mwSize *idims;
    int nsubs=0;
    
    /* Number of pixels */
    int npixels=1;
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("Six inputs are required.");
    } else if(nlhs!=12) {
        mexErrMsgTxt("Twelve outputs are required");
    }
    
    for(i=0; i<6; i++) {
        if(!mxIsSingle(prhs[i])){ mexErrMsgTxt("Inputs must be single"); }
    }
    
    /*  Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    /* Get the sizes of the inputs */
    idims = mxGetDimensions(prhs[0]);
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
    
    /* Assign pointers to each input. */
    Dxx = (float *)mxGetPr(prhs[0]);
    Dxy = (float *)mxGetPr(prhs[1]);
    Dxz = (float *)mxGetPr(prhs[2]);
    Dyy = (float *)mxGetPr(prhs[3]);
    Dyz = (float *)mxGetPr(prhs[4]);
    Dzz = (float *)mxGetPr(prhs[5]);

	/* Assign pointers to each output. */
    plhs[0] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    Deiga = (float *)mxGetPr(plhs[0]); Deigb = (float *)mxGetPr(plhs[1]); Deigc = (float *)mxGetPr(plhs[2]);
    
	/* eigenvectors) */
    plhs[3] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[5] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
	Dvecxa = (float *)mxGetPr(plhs[3]); Dvecya = (float *)mxGetPr(plhs[4]); Dvecza = (float *)mxGetPr(plhs[5]);
    
    plhs[6] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[7] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[8] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
	Dvecxb = (float *)mxGetPr(plhs[6]); Dvecyb = (float *)mxGetPr(plhs[7]); Dveczb = (float *)mxGetPr(plhs[8]);
    
    plhs[9] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[10] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[11] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
	Dvecxc = (float *)mxGetPr(plhs[9]); Dvecyc = (float *)mxGetPr(plhs[10]); Dveczc = (float *)mxGetPr(plhs[11]);
        
	for(i=0; i<npixels; i++) {
		Ma[0][0]=(double)Dxx[i]; Ma[0][1]=(double)Dxy[i]; Ma[0][2]=(double)Dxz[i];
		Ma[1][0]=(double)Dxy[i]; Ma[1][1]=(double)Dyy[i]; Ma[1][2]=(double)Dyz[i];
		Ma[2][0]=(double)Dxz[i]; Ma[2][1]=(double)Dyz[i]; Ma[2][2]=(double)Dzz[i];
		eigen_decomposition(Ma, Davec, Daeig);
        
		Deiga[i]=(float)Daeig[2]; 
        Deigb[i]=(float)Daeig[1]; 
        Deigc[i]=(float)Daeig[0];
        
		Dvecxa[i]=(float)Davec[0][2]; Dvecya[i]=(float)Davec[1][2]; Dvecza[i]=(float)Davec[2][2];
		Dvecxb[i]=(float)Davec[0][1]; Dvecyb[i]=(float)Davec[1][1]; Dveczb[i]=(float)Davec[2][1];
		Dvecxc[i]=(float)Davec[0][0]; Dvecyc[i]=(float)Davec[1][0]; Dveczc[i]=(float)Davec[2][0];
	}
}


