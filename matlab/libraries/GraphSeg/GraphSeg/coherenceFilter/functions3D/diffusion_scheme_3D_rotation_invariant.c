#include "mex.h"
#include "math.h"
#ifndef min
#define min(a,b)        ((a) < (b) ? (a): (b))
#endif
#ifndef max
#define max(a,b)        ((a) > (b) ? (a): (b))
#endif

#ifndef mwSize
#define mwSize int
#endif
__inline int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

void gradient3Dx_float(float *I, int *sizeI, float *Ix)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    float *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const float smoothfilter[3]={0.187500f,0.625000f,0.187500f};
    const float derivafilter[3]={-0.5f,0.0f,0.5f};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(float *)malloc(4*nSlice*sizeof(float));
    Irow=(float *)malloc(sizeI[0]*sizeof(float));
    
    for(z=0; z<sizeI[2]; z++)
    {
        offsetz=nSlice*z;
        offset_slice=nSlice*slice_select;
        
        for(y=0; y<sizeI[1]; y++)
        {
          /* Smooth y - direction  */
            yn=max(y-1,0);
            yp=min(y+1,sizeI[1]-1);
            
            indexn=yn*sizeI[0]+offsetz; 
            indexc=y*sizeI[0]+offsetz;
            indexp=yp*sizeI[0]+offsetz;
            
            for(x=0; x<sizeI[0]; x++)
            {
                Irow[x] =smoothfilter[0]*I[indexn+x];
                Irow[x]+=smoothfilter[1]*I[indexc+x];
                Irow[x]+=smoothfilter[2]*I[indexp+x];
            }

            indexc=y*sizeI[0]+offset_slice;
          /*  Gradient in x - direction  */
            for(x=0; x<sizeI[0]; x++)
            {
                xn=max(x-1,0); xp=min(x+1,sizeI[0]-1);
                Islices[indexc+x]=derivafilter[0]*Irow[xn]+derivafilter[1]*Irow[x]+derivafilter[2]*Irow[xp];
            }
        }
       
        /* Smooth in z - direction  */
        if(z==1) /* Forward          */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            for(i=0; i<nSlice; i++) 
            { 
                Ix[i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
        else if(z>1) /* Central  */
        {
            indexn=slice_select_p2*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*(z-1);
            for(i=0; i<nSlice; i++) 
            { 
                Ix[offsetz+i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
        
        if(z==(sizeI[1]-1)) /* Backward  */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*z;
            for(i=0; i<nSlice; i++) 
            { 
                Ix[offsetz+i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
       
        slice_select_p2=slice_select_p1; slice_select_p1=slice_select; slice_select++; if(slice_select>3) { slice_select=0; }

    }
    free(Irow);
    free(Islices);
}

void gradient3Dy_float(float *I, int *sizeI, float *Iy)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    float *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const float smoothfilter[3]={0.187500f,0.625000f,0.187500f};
    const float derivafilter[3]={-0.5f,0.0f,0.5f};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(float *)malloc(4*nSlice*sizeof(float));
    Irow=(float *)malloc(sizeI[0]*sizeof(float));
    
    for(z=0; z<sizeI[2]; z++)
    {
        offsetz=nSlice*z;
        offset_slice=nSlice*slice_select;
        
        for(y=0; y<sizeI[1]; y++)
        {
          /* Smooth y - direction  */
            yn=max(y-1,0);
            yp=min(y+1,sizeI[1]-1);
            
            indexn=yn*sizeI[0]+offsetz; 
            indexc=y*sizeI[0]+offsetz;
            indexp=yp*sizeI[0]+offsetz;
            
            for(x=0; x<sizeI[0]; x++)
            {
                Irow[x] =derivafilter[0]*I[indexn+x];
                Irow[x]+=derivafilter[1]*I[indexc+x];
                Irow[x]+=derivafilter[2]*I[indexp+x];
            }

            indexc=y*sizeI[0]+offset_slice;
          /*  Gradient in x - direction  */
            for(x=0; x<sizeI[0]; x++)
            {
                xn=max(x-1,0); xp=min(x+1,sizeI[0]-1);
                Islices[indexc+x]=smoothfilter[0]*Irow[xn]+smoothfilter[1]*Irow[x]+smoothfilter[2]*Irow[xp];
            }
        }
       
        /* Smooth in z - direction  */
        if(z==1) /* Forward          */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            for(i=0; i<nSlice; i++) 
            { 
                Iy[i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
        else if(z>1) /* Central  */
        {
            indexn=slice_select_p2*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*(z-1);
            for(i=0; i<nSlice; i++) 
            { 
                Iy[offsetz+i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
        
        if(z==(sizeI[1]-1)) /* Backward  */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*z;
            for(i=0; i<nSlice; i++) 
            { 
                Iy[offsetz+i]=smoothfilter[0]*Islices[i+indexn]+smoothfilter[1]*Islices[i+indexc]+smoothfilter[2]*Islices[i+indexp];
            }
        }
       
        slice_select_p2=slice_select_p1; slice_select_p1=slice_select; slice_select++; if(slice_select>3) { slice_select=0; }
    }
    free(Irow);
    free(Islices);
}

void gradient3Dz_float(float *I, int *sizeI, float *Iz)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    float *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const float smoothfilter[3]={0.187500f,0.625000f,0.187500f};
    const float derivafilter[3]={-0.5f,0.0f,0.5f};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(float *)malloc(4*nSlice*sizeof(float));
    Irow=(float *)malloc(sizeI[0]*sizeof(float));
    
    for(z=0; z<sizeI[2]; z++)
    {
        offsetz=nSlice*z;
        offset_slice=nSlice*slice_select;
        
        for(y=0; y<sizeI[1]; y++)
        {
          /* Smooth y - direction  */
            yn=max(y-1,0);
            yp=min(y+1,sizeI[1]-1);
            
            indexn=yn*sizeI[0]+offsetz; 
            indexc=y*sizeI[0]+offsetz;
            indexp=yp*sizeI[0]+offsetz;
            
            for(x=0; x<sizeI[0]; x++)
            {
                Irow[x] =smoothfilter[0]*I[indexn+x];
                Irow[x]+=smoothfilter[1]*I[indexc+x];
                Irow[x]+=smoothfilter[2]*I[indexp+x];
            }

            indexc=y*sizeI[0]+offset_slice;
          /*  Gradient in x - direction  */
            for(x=0; x<sizeI[0]; x++)
            {
                xn=max(x-1,0); xp=min(x+1,sizeI[0]-1);
                Islices[indexc+x]=smoothfilter[0]*Irow[xn]+smoothfilter[1]*Irow[x]+smoothfilter[2]*Irow[xp];
            }
        }
       
        /* Smooth in z - direction  */
        if(z==1) /* Forward          */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            for(i=0; i<nSlice; i++) 
            { 
                Iz[i]=derivafilter[0]*Islices[i+indexn]+derivafilter[1]*Islices[i+indexc]+derivafilter[2]*Islices[i+indexp];
            }
        }
        else if(z>1) /* Central  */
        {
            indexn=slice_select_p2*nSlice; indexc=slice_select_p1*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*(z-1);
            for(i=0; i<nSlice; i++) 
            { 
                Iz[offsetz+i]=derivafilter[0]*Islices[i+indexn]+derivafilter[1]*Islices[i+indexc]+derivafilter[2]*Islices[i+indexp];
            }
        }
        
        if(z==(sizeI[1]-1)) /* Backward  */
        {
            indexn=slice_select_p1*nSlice; indexc=slice_select*nSlice; indexp=slice_select*nSlice;
            offsetz=nSlice*z;
            for(i=0; i<nSlice; i++) 
            { 
                Iz[offsetz+i]=derivafilter[0]*Islices[i+indexn]+derivafilter[1]*Islices[i+indexc]+derivafilter[2]*Islices[i+indexp];
            }
        }
       
        slice_select_p2=slice_select_p1; slice_select_p1=slice_select; slice_select++; if(slice_select>3) { slice_select=0; }

    }
    free(Irow);
    free(Islices);
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    float *u, *Dxx, *Dyy, *Dzz, *Dxy, *Dxz, *Dyz;
    float *u_new;
	float *j1, *j2, *j3, *ud, *du;
    const mwSize *dimsu_const;
    int dimsu[3], npixels;
    float dt=0.5f;
	double *dt_d;
	int i;
    int x, y,z,index;
   
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
    npixels=dimsu[0]*dimsu[1]*dimsu[2];
	
    /* Assign pointers to each input. */
    u=(float *)mxGetData(prhs[0]); 
    Dxx=(float *)mxGetData(prhs[1]); /* Dxx */
    Dyy=(float *)mxGetData(prhs[4]); /* Dyy */
    Dzz=(float *)mxGetData(prhs[6]); /* Dzz */
    Dxy=(float *)mxGetData(prhs[2]); /* Dxy */
    Dxz=(float *)mxGetData(prhs[3]); /* Dxz */
    Dyz=(float *)mxGetData(prhs[5]); /* Dyz */

 	j1 =(float *)malloc(dimsu[0]*dimsu[1]*dimsu[2]*sizeof(float));
	j2 =(float *)malloc(dimsu[0]*dimsu[1]*dimsu[2]*sizeof(float));
	j3 =(float *)malloc(dimsu[0]*dimsu[1]*dimsu[2]*sizeof(float));
	ud =(float *)malloc(dimsu[0]*dimsu[1]*dimsu[2]*sizeof(float));
	

	/* 3 : Calculate the flux components */
	/* j1 = Dxx .* ux + Dxy .*uy + Dxz .*uz; */
	/* j2 = Dxy .* ux + Dyy .*uy + Dyz .*uz; */
	/* j3 = Dxz .* ux + Dyz .*uy + Dzz .*uz; */

	gradient3Dx_float(u, dimsu, ud);
	for (i=0; i<npixels; i++) 
	{ 
		j1[i]=Dxx[i]*ud[i]; j2[i]=Dxy[i]*ud[i]; j3[i]=Dxz[i]*ud[i];
	}
	gradient3Dy_float(u, dimsu, ud);
	for (i=0; i<npixels; i++) 
	{ 
		j1[i]+=Dxy[i]*ud[i]; j2[i]+=Dyy[i]*ud[i]; j3[i]+=Dyz[i]*ud[i];
	}
	gradient3Dz_float(u, dimsu, ud);
	for (i=0; i<npixels; i++) 
	{ 
		j1[i]+=Dxz[i]*ud[i]; j2[i]+=Dyz[i]*ud[i]; j3[i]+=Dzz[i]*ud[i];
	}
    
    /*j1(:,:,1)=0; j1(:,:,end)=0; j1(:,1,:)=0; j1(:,end,:)=0; j1(1,:,:)=0; j1(end,:,:)=0; */
    /*j2(:,:,1)=0; j2(:,:,end)=0; j2(:,1,:)=0; j2(:,end,:)=0; j2(1,:,:)=0; j2(end,:,:)=0; */
    /*j3(:,:,1)=0; j3(:,:,end)=0; j3(:,1,:)=0; j3(:,end,:)=0; j3(1,:,:)=0; j3(end,:,:)=0; */

    
	for (y=0; y<dimsu[1]; y++) { 
        for (x=0; x<dimsu[0]; x++) { 
            index=mindex3(x,y,0,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0; 
            index=mindex3(x,y,dimsu[2]-1,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0;
        }
    }
        
    for (z=0; z<dimsu[2]; z++) { 
        for (x=0; x<dimsu[0]; x++) { 
            index=mindex3(x,0,z,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0; 
            index=mindex3(x,dimsu[1]-1,z,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0;
        }
    }
    
    for (z=0; z<dimsu[2]; z++) { 
        for (y=0; y<dimsu[1]; y++) { 
            index=mindex3(0,y,z,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0; 
            index=mindex3(dimsu[0]-1,y,z,dimsu[0],dimsu[1]);
            j1[index]=0; j2[index]=0; j3[index]=0;
        }
    }
    
    
    /* 4 : Calculate ... by means of the optimized derivative filters */
	/* du = derivatives(j1,'x')+derivatives(j2,'y')+derivatives(j3,'z'); */
	du =(float *)malloc(dimsu[0]*dimsu[1]*dimsu[2]*sizeof(float));	
	gradient3Dx_float(j1, dimsu, du);
	gradient3Dy_float(j2, dimsu, ud);
	for (i=0; i<npixels; i++) { du[i]+=ud[i]; }
	gradient3Dz_float(j3, dimsu, ud);
	for (i=0; i<npixels; i++) { du[i]+=ud[i]; }
	
	/* Free memory */
	free(j1);
	free(j2);
	free(j3);
	free(ud);
	
	/* 5 : Update in an explicit way. */
	/* u=u+du*dt; */

	/* Create output array */
    plhs[0] = mxCreateNumericArray(3, dimsu, mxSINGLE_CLASS, mxREAL);

	/* Assign pointer to output. */
    u_new = (float *)mxGetData(plhs[0]);
	
	for (i=0; i<npixels; i++) 
	{ 
		u_new[i]+=u[i]+du[i]*dt; 
	}
	
	/* Free memory */
	free(du);
}


