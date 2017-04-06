#include "mex.h"
#include "math.h"
#include "stdlib.h"
#ifndef min
#define min(a,b)        ((a) < (b) ? (a): (b))
#endif
#ifndef max
#define max(a,b)        ((a) > (b) ? (a): (b))
#endif
#ifndef mwSize
#define mwSize int
#endif
#include "CoherenceFilterStep3D_functions.c"

struct options {
    double T;
    double dt;
    double sigma;
    double rho;
    double C;
    double m;
    double alpha;
};

void setdefaultoptions(struct options* t) {
    t->T=2;
    t->dt=0.1;
    t->sigma=1;
    t->rho=1;
    t->C=1e-10;
    t->m=1;
    t->alpha=0.001;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* Input image volume and ouput image volume */
    float *u, *u_new;
    
    /* Options structure variables */
    mxArray *TempField;
    double *OptionsField;
    int field_num;
    struct options Options;
    
    /* Size input image volume */
    int ndimsu;
    const mwSize *dimsu_const;
    int dimsu[3];
    int npixelsu;
    
    float *usigma; /* Gaussian filtered image volume */
    float *ux, *uy, *uz; /* Gradients of smoothed image */
    float *Jxx, *Jxy, *Jxz, *Jyy, *Jyz, *Jzz, *J[6]; /* Structure tensor */
    float *Dxx, *Dxy, *Dxz, *Dyy, *Dyz, *Dzz; /* Diffusion tensor */
    
    
    /* Check number of inputs/outputs */
    if(nrhs<1) { mexErrMsgTxt("1 input variable is required, the other one is optional."); }
    if(nlhs<1) { mexErrMsgTxt("1 output variable is required"); }
        
    /* Check type of inputs */
    if(!mxIsSingle(prhs[0])) { mexErrMsgTxt("Input Image volume must be of datatype Single");}
    if(nrhs==2){ if(!mxIsStruct(prhs[1])){ mexErrMsgTxt("Options must be of type structure"); } }
    
    /* Set Options struct */
    setdefaultoptions(&Options);
    if(nrhs==2) {
        field_num = mxGetFieldNumber(prhs[1], "T");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("aValues in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.T=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[1], "dt");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.dt=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[1], "sigma");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.sigma=OptionsField[0];
        }     
        field_num = mxGetFieldNumber(prhs[1], "rho");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.rho=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[1], "C");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.C=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[1], "m");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.m=OptionsField[0];
        }  
        field_num = mxGetFieldNumber(prhs[1], "alpha");
        if(field_num>=0) {
            TempField=mxGetFieldByNumber(prhs[1], 0, field_num);
            if(!mxIsDouble(TempField)) { mexErrMsgTxt("Values in options structure must be of datatype double"); }
            OptionsField=mxGetPr(TempField);
            Options.alpha=OptionsField[0];
        }  
    }
    
    /* Check and get input image dimensions */
    ndimsu=mxGetNumberOfDimensions(prhs[0]);
    if(ndimsu!=3) { mexErrMsgTxt("Input Image must be 3D"); }
    dimsu_const = mxGetDimensions(prhs[0]);
    dimsu[0]=dimsu_const[0]; dimsu[1]=dimsu_const[1]; dimsu[2]=dimsu_const[2];
    npixelsu=dimsu[0]*dimsu[1]*dimsu[2];
    
    /* Connect input */
    u =(float *)mxGetData(prhs[0]);

    /* Gaussian Filtering of input image volume*/
    usigma = mallocf(npixelsu);  
    GaussianFiltering3D_float(u, usigma, dimsu, Options.sigma, 4*Options.sigma);
    
    /* Calculate the image gradients of the smoothed image volume */
    ux = mallocf(npixelsu);  
    uy = mallocf(npixelsu);  
    uz = mallocf(npixelsu);   
    
    gradient3Dx_float(usigma, dimsu, ux);
    gradient3Dy_float(usigma, dimsu, uy);
    gradient3Dz_float(usigma, dimsu, uz);
    
    /* remove usigma from memory */
    free(usigma);
    
    /* Compute the 3D structure tensors J of the image */
    StructureTensor3D(ux,uy,uz, J, dimsu, Options.rho);
    Jxx=J[0]; Jyy=J[1]; Jzz=J[2]; Jxy=J[3]; Jxz=J[4]; Jyz=J[5];

    /* remove gradients from memory */
    free(ux); free(uy); free(uz);
  
    /* Structure to Diffusion Tensor Weickert */
    Dxx = mallocf(npixelsu);  
    Dyy = mallocf(npixelsu);  
    Dzz = mallocf(npixelsu);  
    Dxy = mallocf(npixelsu);  
    Dxz = mallocf(npixelsu);  
    Dyz = mallocf(npixelsu);  
       
    StructureTensor2DiffusionTensor(Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,dimsu,Options.C,Options.m,Options.alpha); 
    
    /* remove structure tensor from memory */
    free(Jxx); free(Jyy); free(Jzz); free(Jxy); free(Jxz); free(Jyz);
   
    /* Create output array */
    plhs[0] = mxCreateNumericArray(3, dimsu, mxSINGLE_CLASS, mxREAL);

	/* Assign pointer to output. */
    u_new = (float *)mxGetData(plhs[0]);
	
    /* Perform the image diffusion */
    diffusion_scheme_3D_rotation_invariance(u,u_new,dimsu,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
    
    /* remove diffusion tensor from memory */
    free(Dxx); free(Dyy); free(Dzz); free(Dxy); free(Dxz); free(Dyz);
}
