#include "mex.h"
#include "math.h"
#include "EigenDecomposition3.c"

#ifndef mwSize
#define mwSize int
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    float *Dxx, *Dxy, *Dxz, *Dyy, *Dyz, *Dzz;
    float *Jxx, *Jxy, *Jxz, *Jyy, *Jyz, *Jzz;

    /* Matrices of Eigenvector calculation */
    double Ma[3][3];
    double Davec[3][3];
    double Daeig[3];
    
    /* Eigenvector and eigenvalues as scalars */
    double mu1, mu2, mu3, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
    
    /* Amplitudes of diffustion tensor */
    double lambda1, lambda2, lambda3;
      
    /* Constants used in diffusion tensor calculation */
    double alpha, C, m;

    /* Denominator  */
    double di;
    /* Eps for finite values */
    float eps=(float)1e-20;
        

    /* Used to read Options structure */
    double *OptionsField;
    int field_num;
    
    /* Loop variable */
    int i;
    
    /* Size of input */
    const mwSize *idims;
    int nsubs=0;
    
    /* Number of pixels */
    int npixels=1;
    
    /* Check for proper number of arguments. */
    if(nrhs!=7) {
        mexErrMsgTxt("Seven inputs are required.");
    } else if(nlhs!=6) {
        mexErrMsgTxt("Six outputs are required");
    }
    
    for(i=0; i<6; i++) {
        if(!mxIsSingle(prhs[i])){ mexErrMsgTxt("Inputs must be single"); }
    }
    
    /* Get constants from Options structure */
    if(!mxIsStruct(prhs[6])){ mexErrMsgTxt("Options must be structure"); }
    field_num = mxGetFieldNumber(prhs[6], "alpha");
    if(field_num>=0) {
        OptionsField=mxGetPr(mxGetFieldByNumber(prhs[6], 0, field_num));
        alpha=(double)OptionsField[0];
    }
    field_num = mxGetFieldNumber(prhs[6], "C");
    if(field_num>=0) {
        OptionsField=mxGetPr(mxGetFieldByNumber(prhs[6], 0, field_num));
        C=(double)OptionsField[0];
    }
    field_num = mxGetFieldNumber(prhs[6], "m");
    if(field_num>=0) {
        OptionsField=mxGetPr(mxGetFieldByNumber(prhs[6], 0, field_num));
        m=(double)OptionsField[0];
    }
    
    
    /*  Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    /* Get the sizes of the inputs */
    idims = mxGetDimensions(prhs[0]);
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
    
    /* Create output Tensor Volumes */
    plhs[0] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);
    plhs[5] = mxCreateNumericArray(nsubs, idims, mxSINGLE_CLASS, mxREAL);

    /* Assign pointers to each input. */
    Jxx = (float *)mxGetPr(prhs[0]);
    Jxy = (float *)mxGetPr(prhs[1]);
    Jxz = (float *)mxGetPr(prhs[2]);
    Jyy = (float *)mxGetPr(prhs[3]);
    Jyz = (float *)mxGetPr(prhs[4]);
    Jzz = (float *)mxGetPr(prhs[5]);

    
    /* Assign pointers to each output. */
    Dxx = (float *)mxGetPr(plhs[0]); 
    Dxy = (float *)mxGetPr(plhs[1]); 
    Dxz = (float *)mxGetPr(plhs[2]);
    Dyy = (float *)mxGetPr(plhs[3]);
    Dyz = (float *)mxGetPr(plhs[4]);
    Dzz = (float *)mxGetPr(plhs[5]);
        
	for(i=0; i<npixels; i++) {
        /* Calculate eigenvectors and values of local Hessian */
		Ma[0][0]=(double)Jxx[i]+eps; Ma[0][1]=(double)Jxy[i]; Ma[0][2]=(double)Jxz[i];
		Ma[1][0]=(double)Jxy[i]; Ma[1][1]=(double)Jyy[i]+eps; Ma[1][2]=(double)Jyz[i];
		Ma[2][0]=(double)Jxz[i]; Ma[2][1]=(double)Jyz[i]; Ma[2][2]=(double)Jzz[i]+eps;
		eigen_decomposition(Ma, Davec, Daeig);

        /* Convert eigenvector and eigenvalue matrices back to scalar variables */
        mu1=Daeig[2]; 
        mu2=Daeig[1]; 
        mu3=Daeig[0];
		v1x=Davec[0][0]; v1y=Davec[1][0]; v1z=Davec[2][0];
		v2x=Davec[0][1]; v2y=Davec[1][1]; v2z=Davec[2][1];
        v3x=Davec[0][2]; v3y=Davec[1][2]; v3z=Davec[2][2];
        
        /* Scaling of diffusion tensor */
        di=mu1-mu3;
        if((di<eps)&&(di>-eps)) { lambda1 = alpha; } else { lambda1 = alpha + (1.0- alpha)*exp(-C/pow(di,(2.0*m))); }
        lambda2 = alpha;
		lambda3 = alpha;
      
        /* Construct the diffusion tensor */
        Dxx[i] = (float)(lambda1*v1x*v1x + lambda2*v2x*v2x + lambda3*v3x*v3x);
        Dyy[i] = (float)(lambda1*v1y*v1y + lambda2*v2y*v2y + lambda3*v3y*v3y);
        Dzz[i] = (float)(lambda1*v1z*v1z + lambda2*v2z*v2z + lambda3*v3z*v3z);
        Dxy[i] = (float)(lambda1*v1x*v1y + lambda2*v2x*v2y + lambda3*v3x*v3y);
        Dxz[i] = (float)(lambda1*v1x*v1z + lambda2*v2x*v2z + lambda3*v3x*v3z);
        Dyz[i] = (float)(lambda1*v1y*v1z + lambda2*v2y*v2z + lambda3*v3y*v3z);
	}
}






