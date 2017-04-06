#include "mex.h"
#include "math.h"

#ifndef mwSize
#define mwSize int
#endif

/* function u=DoTensorDrivenDiffusion3D(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options)

*/

__inline int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

float CalculateDiffusionNewGreyValue(float *u,float *a,float *b,float *c,float *d,float *e,float *f,int x,int y,int z,int nx,int ny,int nz,int px,int py,int pz, float dt, int *dimsu)
{
	float A2, A4, A5, A6, A8,B1, B2,B3,B4,B5,B6,B7,B8,B9,C2,C4,C5,C6,C8;
 	float u_new;
    float du;
    float di;
    float eps=(float)1e-35;
	int index;

	/* Compute tensor-driven diffusion (as in [1] pp. 80-82) */
	A2 = -f[mindex3(x,y,nz,dimsu[0],dimsu[1])]-f[mindex3(x,py,z,dimsu[0],dimsu[1])]; 
	A4 =  e[mindex3(x,y,nz,dimsu[0],dimsu[1])]+e[mindex3(nx,y,z,dimsu[0],dimsu[1])];
	A6 = -e[mindex3(x,y,nz,dimsu[0],dimsu[1])]-e[mindex3(px,y,z,dimsu[0],dimsu[1])];
	A8 =  f[mindex3(x,y,nz,dimsu[0],dimsu[1])]+f[mindex3(x,ny,z,dimsu[0],dimsu[1])];
	B1 = -d[mindex3(nx,y,z,dimsu[0],dimsu[1])]-d[mindex3(x,py,z,dimsu[0],dimsu[1])];
	B3 =  d[mindex3(px,y,z,dimsu[0],dimsu[1])]+d[mindex3(x,py,z,dimsu[0],dimsu[1])];
	B7 =  d[mindex3(nx,y,z,dimsu[0],dimsu[1])]+d[mindex3(x,ny,z,dimsu[0],dimsu[1])];
	B9 = -d[mindex3(px,y,z,dimsu[0],dimsu[1])]-d[mindex3(x,ny,z,dimsu[0],dimsu[1])];
	C2 =  f[mindex3(x,y,pz,dimsu[0],dimsu[1])]+f[mindex3(x,py,z,dimsu[0],dimsu[1])];
	C4 = -e[mindex3(x,y,pz,dimsu[0],dimsu[1])]-e[mindex3(nx,y,z,dimsu[0],dimsu[1])];
	C6 =  e[mindex3(x,y,pz,dimsu[0],dimsu[1])]+e[mindex3(px,y,z,dimsu[0],dimsu[1])];
	C8 = -f[mindex3(x,y,pz,dimsu[0],dimsu[1])]-f[mindex3(x,ny,z,dimsu[0],dimsu[1])];
	
	A5 = c[mindex3(x,y,nz,dimsu[0],dimsu[1])]+c[mindex3(x,y,z,dimsu[0],dimsu[1])];
	B2 = b[mindex3(x,py,z,dimsu[0],dimsu[1])]+b[mindex3(x,y,z,dimsu[0],dimsu[1])];
	B4 = a[mindex3(nx,y,z,dimsu[0],dimsu[1])]+a[mindex3(x,y,z,dimsu[0],dimsu[1])];
	B5 =  -(a[mindex3(nx,y,z,dimsu[0],dimsu[1])] + 2*a[mindex3(x,y,z,dimsu[0],dimsu[1])] + a[mindex3(px,y,z,dimsu[0],dimsu[1])]);
	B5 += -(b[mindex3(x,ny,z,dimsu[0],dimsu[1])] + 2*b[mindex3(x,y,z,dimsu[0],dimsu[1])] + b[mindex3(x,py,z,dimsu[0],dimsu[1])]);
	B5 += -(c[mindex3(x,y,nz,dimsu[0],dimsu[1])] + 2*c[mindex3(x,y,z,dimsu[0],dimsu[1])] + c[mindex3(x,y,pz,dimsu[0],dimsu[1])]);
	B6 = a[mindex3(px,y,z,dimsu[0],dimsu[1])]+a[mindex3(x,y,z,dimsu[0],dimsu[1])];
	B8 = b[mindex3(x,ny,z,dimsu[0],dimsu[1])]+b[mindex3(x,y,z,dimsu[0],dimsu[1])];
	C5 = c[mindex3(x,y,pz,dimsu[0],dimsu[1])]+c[mindex3(x,y,z,dimsu[0],dimsu[1])];

	du=0;
	du+=A2*u[mindex3(x,py,nz,dimsu[0],dimsu[1])];
	du+=A4*u[mindex3(nx,y,nz,dimsu[0],dimsu[1])]; 
	du+=A6*u[mindex3(px,y,nz,dimsu[0],dimsu[1])]; 
	du+=A8*u[mindex3(x,ny,nz,dimsu[0],dimsu[1])]; 
	
	du+=B1*u[mindex3(nx,py,z,dimsu[0],dimsu[1])]; 
	du+=B3*u[mindex3(px,py,z,dimsu[0],dimsu[1])]; 
	du+=B7*u[mindex3(nx,ny,z,dimsu[0],dimsu[1])]; 
	du+=B9*u[mindex3(px,ny,z,dimsu[0],dimsu[1])]; 
	
	du+=C2*u[mindex3(x,py,pz,dimsu[0],dimsu[1])];
	du+=C4*u[mindex3(nx,y,pz,dimsu[0],dimsu[1])]; 
	du+=C6*u[mindex3(px,y,pz,dimsu[0],dimsu[1])]; 
	du+=C8*u[mindex3(x,ny,pz,dimsu[0],dimsu[1])];

	du*=0.5;
	
	du+=A5*u[mindex3(x,y,nz,dimsu[0],dimsu[1])];  
	du+=B2*u[mindex3(x,py,z,dimsu[0],dimsu[1])]; 
	du+=B4*u[mindex3(nx,y,z,dimsu[0],dimsu[1])]; 
	du+=B6*u[mindex3(px,y,z,dimsu[0],dimsu[1])]; 
	du+=B8*u[mindex3(x,ny,z,dimsu[0],dimsu[1])];  
	du+=C5*u[mindex3(x,y,pz,dimsu[0],dimsu[1])];  
	
	du*=0.5f*dt;

	/* Perform the edge preserving diffusion filtering on the image */
	index = mindex3(x,y,z,dimsu[0],dimsu[1]);
	di=(1-0.5f*dt*B5); if((di<eps)&&(di>-eps)) { di=eps; }
	u_new = (u[index] + du )/di;  
	return u_new;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    float *u, *a, *b, *c, *d, *e, *f;
    float *u_new;
    const mwSize *dimsu_const;
    int dimsu[3];
    float dt=0.5f;
	double *dt_d;
    int index;
    int i;
    int x,y,z, nx, ny, nz, px, py, pz;
    
    /* Check number of inputs */
    if(nrhs<8) { mexErrMsgTxt("8 input variables are required."); }
    
    /* Check if inputs are of Single type */
    for(i=0; i<7; i++) {
        if(!mxIsSingle(prhs[i])){ mexErrMsgTxt("Inputs u,Dxx,Dxy,Dxz,Dyy,Dyz and Dzz must be single"); }
    }    
    
    /* Get time stepsize */
	if(!mxIsDouble(prhs[7])){ mexErrMsgTxt("Input dt must be double"); }
    dt_d=(double *)mxGetData(prhs[7]); 
	dt=(float) dt_d[0];
	
    
     /* Check properties of image I */
    if(mxGetNumberOfDimensions(prhs[0])!=3) { mexErrMsgTxt("Image must be 3D"); }
    /* Get the sizes of the image */
    dimsu_const = mxGetDimensions(prhs[0]);
	dimsu[0]=dimsu_const[0];
    dimsu[1]=dimsu_const[1];
    dimsu[2]=dimsu_const[2];
    
    /* Assign pointers to each input. */
    u=(float *)mxGetData(prhs[0]); 
    a=(float *)mxGetData(prhs[1]); /* Dxx */
    b=(float *)mxGetData(prhs[4]); /* Dyy */
    c=(float *)mxGetData(prhs[6]); /* Dzz */
    d=(float *)mxGetData(prhs[2]); /* Dxy */
    e=(float *)mxGetData(prhs[3]); /* Dxz */
    f=(float *)mxGetData(prhs[5]); /* Dyz */

    /* Create output array */
    plhs[0] = mxCreateNumericArray(3, dimsu, mxSINGLE_CLASS, mxREAL);
    /* Assign pointer to output. */
    u_new = (float *)mxGetData(plhs[0]);

    
    
    /* Loop through all pixels in the volume */
    for (z=0; z<dimsu[2]; z++) {
        /* Neighbor coordinates */
        nz=z-1; if(nz<0) {nz=0; }
        pz=z+1; if(pz>(dimsu[2]-1)) { pz=dimsu[2]-1; }
        for (y=0; y<dimsu[1]; y++) {
            ny=y-1; if(ny<0) {ny=0; }
            py=y+1; if(py>(dimsu[1]-1)) { py=dimsu[1]-1; }
            for (x=0; x<dimsu[0]; x++) {
                nx=x-1; if(nx<0) {nx=0; }
                px=x+1; if(px>(dimsu[0]-1)) { px=dimsu[0]-1; }

				index = mindex3(x,y,z,dimsu[0],dimsu[1]);
				u_new[index] = CalculateDiffusionNewGreyValue(u,a,b,c,d,e,f,x,y,z,nx,ny,nz,px,py,pz,dt,dimsu);
            }
        }
    }

    
    
    
}
                
        
        
