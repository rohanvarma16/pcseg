#include "mex.h"
#include "math.h"
#include "stdlib.h"
#ifndef min
#define min(a,b)        ((a) < (b) ? (a): (b))
#endif
#ifndef max
#define max(a,b)        ((a) > (b) ? (a): (b))
#endif
#include "CoherenceFilterStep2D_functions.c"
#define mwSize int

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
    /* Input image  and ouput image  */
    double *u, *u_new;
    
    /* Options structure variables */
    mxArray *TempField;
    double *OptionsField;
    int field_num;
    struct options Options;
    
    /* Size input image  */
    int ndimsu;
    const mwSize *dimsu_const;
    int dimsu[3];
    int npixels2;
	int npixels3;
	
    
    double *usigma; /* Gaussian filtered image  */
    double *ux, *uy; /* Gradients of smoothed image */
    double *Jxx, *Jxy, *Jyy, *J[3]; /* Structure tensor */
    double *Dxx, *Dxy, *Dyy; /* Diffusion tensor */
    
    
    /* Check number of inputs/outputs */
    if(nrhs<1) { mexErrMsgTxt("1 input variable is required, the other one is optional."); }
    if(nlhs<1) { mexErrMsgTxt("1 output variable is required"); }
        
    /* Check type of inputs */
    if(!mxIsDouble(prhs[0])) { mexErrMsgTxt("Input Image  must be of datatype Double");}
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
	if((ndimsu<2)||(ndimsu>3)) { mexErrMsgTxt("Input Image must be 2D"); }
	
    dimsu_const = mxGetDimensions(prhs[0]);
    dimsu[0]=dimsu_const[0]; dimsu[1]=dimsu_const[1]; 
	if(ndimsu==3) { dimsu[2]=dimsu_const[2]; } else { dimsu[2]=1; }
	npixels2=dimsu[0]*dimsu[1];
    npixels3=dimsu[0]*dimsu[1]*dimsu[2];
    
    /* Connect input */
    u =(double *)mxGetData(prhs[0]);

    /* Gaussian Filtering of input image */
    usigma = mallocd(npixels3);  
	GaussianFiltering2Dcolor_double(u, usigma, dimsu, Options.sigma, 4*Options.sigma);
	    
    /* Calculate the image gradients of the smoothed image  */
    ux = mallocd(npixels3);  
    uy = mallocd(npixels3);  

    gradient2Dx_color(usigma, dimsu, ux);
    gradient2Dy_color(usigma, dimsu, uy);
    
    /* remove usigma from memory */
    free(usigma);
    
    /* Compute the 2D structure tensors J of the image */
    StructureTensor2D(ux,uy, J, dimsu, Options.rho);
    Jxx=J[0]; Jyy=J[1] ;Jxy=J[2];

    /* remove gradients from memory */
    free(ux); free(uy); 
	
    /* Structure to Diffusion Tensor Weickert */
    Dxx = mallocd(npixels2);  
    Dyy = mallocd(npixels2);  
    Dxy = mallocd(npixels2);  
    StructureTensor2DiffusionTensor(Jxx,Jxy,Jyy,Dxx,Dxy,Dyy,dimsu,Options.C,Options.m,Options.alpha); 
    
    /* remove structure tensor from memory */
    free(Jxx); free(Jyy);  free(Jxy); 
	
    /* Create output array */
	if(dimsu[2]==1) { plhs[0] = mxCreateNumericArray(2, dimsu, mxDOUBLE_CLASS, mxREAL); }
	else { plhs[0] = mxCreateNumericArray(3, dimsu, mxDOUBLE_CLASS, mxREAL); }
	
	/* Assign pointer to output. */
    u_new = (double *)mxGetData(plhs[0]);
	
    /* Perform the image diffusion */
    diffusion_scheme_2D_rotation_invariance(u,u_new,dimsu,Dxx,Dxy,Dyy,Options.dt);
    
    /* remove diffusion tensor from memory */
    free(Dxx); free(Dyy); free(Dxy);
}
