#include "mex.h"
#include "math.h"
#define absd(a) ((a)>(-a)?(a):(-a))
#ifndef mwSize
#define mwSize int
#endif
/* function u=DoTensorDrivenDiffusion3D(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options)

*/

__inline int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

float CalculateDiffusionNewGreyValue(float *u,float *Dxx,float *Dyy,float *Dzz,float *Dxy,float *Dxz,float *Dyz,int x,int y,int z,int nx,int ny,int nz,int px,int py,int pz, float dt, int *dimsu)
{
    float WE, WW, WS, WN, WB, WF, WSE, WNW, WNE, WSW, WSF, WNF, WEF, WWF, WSB, WNB, WEB, WWB;
    float Rxx, Ryy, Rzz, Rxy, Rxz, Ryz;
 	float HX=1;
    float HY=1;
    float HZ=1;
    float u_new;
    int index;
    
	int index_nx_ny_nz,	index_x_ny_nz,	index_px_ny_nz,	index_nx_y_nz,	index_x_y_nz,	index_px_y_nz,	index_nx_py_nz,	index_x_py_nz, index_px_py_nz;
	int index_nx_ny_z,	index_x_ny_z,	index_px_ny_z,	index_nx_y_z,	index_px_y_z,	index_nx_py_z,	index_x_py_z, index_px_py_z;
	int index_nx_ny_pz,	index_x_ny_pz,	index_px_ny_pz,	index_nx_y_pz,	index_x_y_pz,	index_px_y_pz,	index_nx_py_pz,	index_x_py_pz, index_px_py_pz;
	
	
	
    index=mindex3(x,y,z,dimsu[0],dimsu[1]);
    
	index_nx_ny_nz=mindex3(nx,ny,nz,dimsu[0],dimsu[1]);
	index_x_ny_nz=mindex3(x,ny,nz,dimsu[0],dimsu[1]);
	index_px_ny_nz=mindex3(px,ny,nz,dimsu[0],dimsu[1]);
	index_nx_y_nz=mindex3(nx,y,nz,dimsu[0],dimsu[1]);
	index_x_y_nz=mindex3(x,y,nz,dimsu[0],dimsu[1]);
	index_px_y_nz=mindex3(px,y,nz,dimsu[0],dimsu[1]);
	index_nx_py_nz=mindex3(nx,py,nz,dimsu[0],dimsu[1]);
	index_x_py_nz=mindex3(x,py,nz,dimsu[0],dimsu[1]);
	index_px_py_nz=mindex3(px,py,nz,dimsu[0],dimsu[1]);
	
	index_nx_ny_z=mindex3(nx,ny,z,dimsu[0],dimsu[1]);
	index_x_ny_z=mindex3(x,ny,z,dimsu[0],dimsu[1]);
	index_px_ny_z=mindex3(px,ny,z,dimsu[0],dimsu[1]);
	index_nx_y_z=mindex3(nx,y,z,dimsu[0],dimsu[1]);
	index_px_y_z=mindex3(px,y,z,dimsu[0],dimsu[1]);
	index_nx_py_z=mindex3(nx,py,z,dimsu[0],dimsu[1]);
	index_x_py_z=mindex3(x,py,z,dimsu[0],dimsu[1]);
	index_px_py_z=mindex3(px,py,z,dimsu[0],dimsu[1]);

	index_nx_ny_pz=mindex3(nx,ny,pz,dimsu[0],dimsu[1]);
	index_x_ny_pz=mindex3(x,ny,pz,dimsu[0],dimsu[1]);
	index_px_ny_pz=mindex3(px,ny,pz,dimsu[0],dimsu[1]);
	index_nx_y_pz=mindex3(nx,y,pz,dimsu[0],dimsu[1]);
	index_x_y_pz=mindex3(x,y,pz,dimsu[0],dimsu[1]);
	index_px_y_pz=mindex3(px,y,pz,dimsu[0],dimsu[1]);
	index_nx_py_pz=mindex3(nx,py,pz,dimsu[0],dimsu[1]);
	index_x_py_pz=mindex3(x,py,pz,dimsu[0],dimsu[1]);
	index_px_py_pz=mindex3(px,py,pz,dimsu[0],dimsu[1]);
	
    Rxx  = dt / (2.0f * HX * HX);
    Ryy  = dt / (2.0f * HY * HY);
    Rzz  = dt / (2.0f * HZ * HZ);
    Rxy  = dt / (4.0f * HX * HY);
    Rxz  = dt / (4.0f * HX * HZ);
    Ryz  = dt / (4.0f * HY * HZ);

    /* WEIGHTS */
    WE  = Rxx * (Dxx[index_px_y_z] + Dxx[index]) - Rxy * (absd(Dxy[index_px_y_z]) + absd(Dxy[index])) - Rxz * (absd(Dxz[index_px_y_z]) + absd(Dxz[index]));
    WW  = Rxx * (Dxx[index_nx_y_z] + Dxx[index]) - Rxy * (absd(Dxy[index_nx_y_z]) + absd(Dxy[index])) - Rxz * (absd(Dxz[index_nx_y_z]) + absd(Dxz[index]));
    WS  = Ryy * (Dyy[index_x_py_z] + Dyy[index]) - Rxy * (absd(Dxy[index_x_py_z]) + absd(Dxy[index])) - Ryz * (absd(Dyz[index_x_py_z]) + absd(Dyz[index]));
    WN  = Ryy * (Dyy[index_x_ny_z] + Dyy[index]) - Rxy * (absd(Dxy[index_x_ny_z]) + absd(Dxy[index])) - Ryz * (absd(Dyz[index_x_ny_z]) + absd(Dyz[index]));
    WB  = Rzz * (Dzz[index_x_y_nz] + Dzz[index]) - Ryz * (absd(Dyz[index_x_y_nz]) + absd(Dyz[index])) - Rxz * (absd(Dxz[index_x_y_nz]) + absd(Dxz[index]));
    WF  = Rzz * (Dzz[index_x_y_pz] + Dzz[index]) - Ryz * (absd(Dyz[index_x_y_pz]) + absd(Dyz[index])) - Rxz * (absd(Dxz[index_x_y_pz]) + absd(Dxz[index]));
    WSE = Rxy * (   Dxy[index_px_py_z] + Dxy[index] + absd(Dxy[index_px_py_z]) + absd(Dxy[index]));
    WNW = Rxy * (   Dxy[index_nx_ny_z] + Dxy[index] + absd(Dxy[index_nx_ny_z]) + absd(Dxy[index]));
    WNE = Rxy * ( - Dxy[index_px_ny_z] - Dxy[index] + absd(Dxy[index_px_ny_z]) + absd(Dxy[index]));
    WSW = Rxy * ( - Dxy[index_nx_py_z] - Dxy[index] + absd(Dxy[index_nx_py_z]) + absd(Dxy[index]));
    WSF = Ryz * (   Dyz[index_x_py_pz] + Dyz[index] + absd(Dyz[index_x_py_pz]) + absd(Dyz[index]));
    WNF = Ryz * ( - Dyz[index_x_ny_pz] - Dyz[index] + absd(Dyz[index_x_ny_pz]) + absd(Dyz[index]));
    WEF = Rxz * (   Dxz[index_px_y_pz] + Dxz[index] + absd(Dxz[index_px_y_pz]) + absd(Dxz[index]));
    WWF = Rxz * ( - Dxz[index_nx_y_pz] - Dxz[index] + absd(Dxz[index_nx_y_pz]) + absd(Dxz[index]));
    WSB = Ryz * ( - Dyz[index_x_py_nz] - Dyz[index] + absd(Dyz[index_x_py_nz]) + absd(Dyz[index]));
    WNB = Ryz * (   Dyz[index_x_ny_nz] + Dyz[index] + absd(Dyz[index_x_ny_nz]) + absd(Dyz[index]));
    WEB = Rxz * ( - Dxz[index_px_y_nz] - Dxz[index] + absd(Dxz[index_px_y_nz]) + absd(Dxz[index]));
    WWB = Rxz * (   Dxz[index_nx_y_nz] + Dxz[index] + absd(Dxz[index_nx_y_nz]) + absd(Dxz[index]));

    /* EVOLUTION */
    u_new = u[index]; 
    u_new += WE  *(u[index_px_y_z]  - u[index]);
    u_new += WW  *(u[index_nx_y_z]  - u[index]);
    u_new += WS  *(u[index_x_py_z]  - u[index]);
	u_new += WN  *(u[index_x_ny_z]  - u[index]);
    u_new += WB  *(u[index_x_y_nz]  - u[index]);
    u_new += WF  *(u[index_x_y_pz]  - u[index]);
	u_new += WSE *(u[index_px_py_z] - u[index]);
    u_new += WNW *(u[index_nx_ny_z] - u[index]);
	u_new += WSW *(u[index_nx_py_z] - u[index]);
    u_new += WNE *(u[index_px_ny_z] - u[index]);
    u_new += WNB *(u[index_x_ny_nz] - u[index]);
    u_new += WNF *(u[index_x_ny_pz] - u[index]);
    u_new += WEB *(u[index_px_y_nz] - u[index]);
    u_new += WEF *(u[index_px_y_pz] - u[index]);
    u_new += WWB *(u[index_nx_y_nz] - u[index]);
    u_new += WWF *(u[index_nx_y_pz] - u[index]);
    u_new += WSB *(u[index_x_py_nz] - u[index]);
    u_new += WSF *(u[index_x_py_pz] - u[index]);
    
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
                
        
        
