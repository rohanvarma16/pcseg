#include "mex.h"
#include "math.h"
/*   undef needed for LCC compiler  */
#undef EXTERN_C
#ifdef _WIN32
	#include <windows.h>
	#include <process.h>
#else
	#include <pthread.h>
#endif
#define clamp(a, b1, b2) min(max(a, b1), b2);
#define absd(a) ((a)>(-a)?(a):(-a))
#define pow2(a) a*a
#define n 3
#define inv3 0.3333333333333333
#define root3 1.7320508075688772

__inline int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}

float * mallocf(int a)
{
    float *ptr = malloc(a*sizeof(float));
    if (ptr == NULL) { mexErrMsgTxt("Out of Memory"); }
    return ptr;
}

void imfilter1D_float(float *I, int lengthI, float *H, int lengthH, float *J) {
    int x, i, index, offset;
    int b2, offset2;
    if(lengthI==1)  
    { 
        J[0]=I[0];
    }
    else
    {
        offset=(lengthH-1)/2;
        for(x=0; x<min(offset,lengthI); x++) {
            J[x]=0;
            b2=lengthI-1; offset2=x-offset;
            for(i=0; i<lengthH; i++) {
                index=clamp(i+offset2, 0, b2); J[x]+=I[index]*H[i];
            }
        }
       
        for(x=offset; x<(lengthI-offset); x++) {
            J[x]=0;
            b2=lengthI-1; offset2=x-offset;
            for(i=0; i<lengthH; i++) {
                index=i+offset2; J[x]+=I[index]*H[i];
            }
        }
       
         b2=lengthI-1; 
         for(x=max(lengthI-offset,offset); x<lengthI; x++) {
              J[x]=0;
              offset2=x-offset;
              for(i=0; i<lengthH; i++) {
                  index=clamp(i+offset2, 0, b2); J[x]+=I[index]*H[i];
             }
         }
       
    }
}

void imfilter2D_float(float *I, int * sizeI, float *H, int lengthH, float *J) {
    int y, x, i, y2;
    float *Irow, *Crow;
    int index=0, line=0;
    float *RCache;
    int *nCache;
    int hks, offset, offset2;
    RCache=mallocf(lengthH*sizeI[0]);
    for(i=0; i<lengthH*sizeI[0]; i++) { RCache[i]=0; }
    nCache=(int *)malloc(lengthH*sizeof(int));
    for(i=0; i<lengthH; i++) { nCache[i]=0; }
    hks=((lengthH-1)/2);
    for(y=0; y<min(hks,sizeI[1]); y++) {
        Irow=&I[index];
        Crow=&RCache[line*sizeI[0]];
        imfilter1D_float(Irow, sizeI[0], H, lengthH, Crow);
        index+=sizeI[0];
        if(y!=(sizeI[1]-1))
        {
            line++; if(line>(lengthH-1)) { line=0; }
        }
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
    for(y2=y; y2<hks; y2++) {
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
            
    for(y=hks; y<(sizeI[1]-1); y++) {
        Irow=&I[index];
        Crow=&RCache[line*sizeI[0]];
        imfilter1D_float(Irow, sizeI[0], H, lengthH, Crow);
        offset=(y-hks)*sizeI[0]; offset2=nCache[0]*sizeI[0];
        for(x=0; x<sizeI[0]; x++) { J[offset+x]=RCache[offset2+x]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*sizeI[0];
            for(x=0; x<sizeI[0]; x++) { J[offset+x]+=RCache[offset2+x]*H[i]; }
        }
        index+=sizeI[0];
        line++; if(line>(lengthH-1)) { line=0; }
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }

    for(y=max(sizeI[1]-1,hks); y<sizeI[1]; y++) {
        Irow=&I[index];
        Crow=&RCache[line*sizeI[0]];
        imfilter1D_float(Irow, sizeI[0], H, lengthH, Crow);
        offset=(y-hks)*sizeI[0]; offset2=nCache[0]*sizeI[0];
        for(x=0; x<sizeI[0]; x++) { J[offset+x]=RCache[offset2+x]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*sizeI[0];
            for(x=0; x<sizeI[0]; x++) { J[offset+x]+=RCache[offset2+x]*H[i]; }
        }
        index+=sizeI[0];
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }

    for(y=max(sizeI[1],hks); y<(sizeI[1]+hks); y++) {
        offset=(y-hks)*sizeI[0]; offset2=nCache[0]*sizeI[0];
        for(x=0; x<sizeI[0]; x++) { J[offset+x]=RCache[offset2+x]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*sizeI[0];
            for(x=0; x<sizeI[0]; x++) { J[offset+x]+=RCache[offset2+x]*H[i]; }
        }
        index+=sizeI[0];
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }

    free(RCache);
}

void imfilter3D_float(float *I, int * sizeI, float *H, int lengthH, float *J) {
    int z, j, i, z2;
    float *Islice, *Cslice;
    int index=0, line=0;
    float *SCache;
    int *nCache;
    int hks, offset, offset2;
    int nslice;
    nslice=sizeI[0]*sizeI[1];
    SCache=mallocf(lengthH*nslice);
	for(i=0; i<nslice; i++) { SCache[i]=0; }
    nCache=(int *)malloc(lengthH*sizeof(int));
    for(i=0; i<lengthH; i++) { nCache[i]=0; }
    hks=((lengthH-1)/2);
    for(z=0; z<min(hks,sizeI[2]); z++) {
        Islice=&I[index];
        Cslice=&SCache[line*nslice];
        imfilter2D_float(Islice, sizeI, H, lengthH, Cslice);
        index+=nslice;
        if(z!=(sizeI[2]-1))
        {
            line++; if(line>(lengthH-1)) { line=0; }
        }
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
    for(z2=z; z2<hks; z2++) {
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
    for(z=hks; z<(sizeI[2]-1); z++) {
        Islice=&I[index];
        Cslice=&SCache[line*nslice];
        imfilter2D_float(Islice, sizeI, H, lengthH, Cslice);
        offset=(z-hks)*nslice; offset2=nCache[0]*nslice;
        for(j=0; j<nslice; j++) { J[offset+j]=SCache[offset2+j]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*nslice;
            for(j=0; j<nslice; j++) { J[offset+j]+=SCache[offset2+j]*H[i]; }
        }
        index+=nslice;
        line++; if(line>(lengthH-1)) { line=0; }
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
    for(z=max(sizeI[2]-1,hks); z<sizeI[2]; z++) {
        Islice=&I[index];
        Cslice=&SCache[line*nslice];
        imfilter2D_float(Islice, sizeI, H, lengthH, Cslice);
        offset=(z-hks)*nslice; offset2=nCache[0]*nslice;
        for(j=0; j<nslice; j++) { J[offset+j]=SCache[offset2+j]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*nslice;
            for(j=0; j<nslice; j++) { J[offset+j]+=SCache[offset2+j]*H[i]; }
        }
        index+=nslice;
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }
    for(z=max(sizeI[2],hks); z<(sizeI[2]+hks); z++) {
        offset=(z-hks)*nslice; offset2=nCache[0]*nslice;
        for(j=0; j<nslice; j++) { J[offset+j]=SCache[offset2+j]*H[0]; }
        for(i=1; i<lengthH; i++) {
            offset2=nCache[i]*nslice;
            for(j=0; j<nslice; j++) { J[offset+j]+=SCache[offset2+j]*H[i]; }
        }
        index+=nslice;
        for(i=0; i<(lengthH-1); i++) { nCache[i]=nCache[i+1]; } nCache[lengthH-1]=line;
    }

    free(SCache);
}

void GaussianFiltering3D_float(float *I, float *J, int *dimsI, double sigma, double kernel_size)
{
	int kernel_length,i;
    double x;
    float *H, totalH=0;
	
	/* Construct the 1D gaussian kernel */
	if(kernel_size<1) { kernel_size=1; }
    kernel_length=(int)(2*ceil(kernel_size/2)+1);
	H = mallocf(kernel_length);
	x=-ceil(kernel_size/2);
	for (i=0; i<kernel_length; i++) { H[i]=(float)exp(-((x*x)/(2*(sigma*sigma)))); totalH+=H[i]; x++; }
	for (i=0; i<kernel_length; i++) { H[i]/=totalH; }
	
	/* Do the filtering */
	imfilter3D_float(I, dimsI, H, kernel_length, J);
    /* Clear memory gaussian kernel */
	free(H);
}


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
    Islices=mallocf(4*nSlice);
    Irow=mallocf(sizeI[0]);
    
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
    Islices=mallocf(4*nSlice);
    Irow=mallocf(sizeI[0]);
    
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
    Islices=mallocf(4*nSlice);
    Irow=mallocf(sizeI[0]);
    
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


void StructureTensor3D(float *ux, float *uy, float *uz,  float **J, int *dimsu, double rho) {
    int npixelsu;
    int i;
    float *Jxx, *Jxy, *Jxz, *Jyy, *Jyz, *Jzz; /* Structure tensor */
       
    npixelsu=dimsu[0]*dimsu[1]*dimsu[2];
    
    Jxx = mallocf(npixelsu);  
    Jyy = mallocf(npixelsu); 
    Jzz = mallocf(npixelsu); 
    Jxy = mallocf(npixelsu);  
    Jxz = mallocf(npixelsu); 
    Jyz = mallocf(npixelsu); 
    
    /* J(grad u_sigma) */
    for(i=0; i<npixelsu; i++)
    {
        Jxx[i] = ux[i]*ux[i]; Jxy[i] = ux[i]*uy[i];
        Jxz[i] = ux[i]*uz[i]; Jyy[i] = uy[i]*uy[i];
        Jyz[i] = uy[i]*uz[i]; Jzz[i] = uz[i]*uz[i];
    }

    /* Do the gaussian smoothing */
    J[0]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jxx, J[0], dimsu, rho, 4*rho);
    free(Jxx);
    J[1]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jyy, J[1], dimsu, rho, 4*rho);
    free(Jyy);
    J[2]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jzz, J[2], dimsu, rho, 4*rho);
    free(Jzz);
    J[3]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jxy, J[3], dimsu, rho, 4*rho);
    free(Jxy);
    J[4]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jxz, J[4], dimsu, rho, 4*rho);
    free(Jxz);
    J[5]=mallocf(npixelsu); 
    GaussianFiltering3D_float(Jyz, J[5], dimsu, rho, 4*rho);
    free(Jyz);
}


/* domain Java Matrix library JAMA. */
static double hypot2(double x, double y) { return sqrt(x*x+y*y); }

/* Symmetric Householder reduction to tridiagonal form. */
static void tred2(double V[n][n], double d[n], double e[n]) {
    
/*  This is derived from the Algol procedures tred2 by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */
    int i, j, k;
    double scale;
    double f, g, h;
    double hh;
    for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }
    
    /* Householder reduction to tridiagonal form. */
    
    for (i = n-1; i > 0; i--) {
        /* Scale to avoid under/overflow. */
        scale = 0.0;
        h = 0.0;
        for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
        } else {
            
            /* Generate Householder vector. */
            
            for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
            f = d[i-1];
            g = sqrt(h);
            if (f > 0) { g = -g; }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++) { e[j] = 0.0; }
            
            /* Apply similarity transformation to remaining columns. */
            
            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
            hh = f / (h + h);
            for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
            for (j = 0; j < i; j++) {
                f = d[j]; g = e[j];
                for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    
    /* Accumulate transformations. */
    
    for (i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
            for (j = 0; j <= i; j++) {
                g = 0.0;
                for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
            }
        }
        for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
    }
    for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
}

/* Symmetric tridiagonal QL algorithm. */
static void tql2(double V[n][n], double d[n], double e[n]) {
    
/*  This is derived from the Algol procedures tql2, by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */
    
    int i, j, k, l, m;
    double f;
    double tst1;
    double eps;
    int iter;
    double g, p, r;
    double dl1, h, c, c2, c3, el1, s, s2;
    
    for (i = 1; i < n; i++) { e[i-1] = e[i]; }
    e[n-1] = 0.0;
    
    f = 0.0;
    tst1 = 0.0;
    eps = pow(2.0, -52.0);
    for (l = 0; l < n; l++) {
        
        /* Find small subdiagonal element */
        
        tst1 = max(tst1, fabs(d[l]) + fabs(e[l]));
        m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) { break; }
            m++;
        }
        
        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */
        
        if (m > l) {
            iter = 0;
            do {
                iter = iter + 1;  /* (Could check iteration count here.) */
                /* Compute implicit shift */
                g = d[l];
                p = (d[l+1] - g) / (2.0 * e[l]);
                r = hypot2(p, 1.0);
                if (p < 0) { r = -r; }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                dl1 = d[l+1];
                h = g - d[l];
                for (i = l+2; i < n; i++) { d[i] -= h; }
                f = f + h;
                /* Implicit QL transformation. */
                p = d[m]; c = 1.0; c2 = c; c3 = c;
                el1 = e[l+1]; s = 0.0; s2 = 0.0;
                for (i = m-1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);
                    /* Accumulate transformation. */
                    for (k = 0; k < n; k++) {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                
                /* Check for convergence. */
            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
    
    /* Sort eigenvalues and corresponding vectors. */
    for (i = 0; i < n-1; i++) {
        k = i;
        p = d[i];
        for (j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}



void roots3(double d[3], double c0,double c1, double c2)
{
    double c2Div3, aDiv3, mbDiv2, q, magnitude, angle, cs, sn;
    
    /* Solve the roots of  y^3 + c2 * y^2 + c1 *y + c0  */
    c2Div3 = -c2*inv3;
    aDiv3 = (c1 + c2*c2Div3)*inv3;
    if (aDiv3 > 0.0) { aDiv3 = 0.0; }
    mbDiv2 = 0.5*(-c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
    q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
    if (q > 0.0) { q = 0.0; }
    magnitude = sqrt(-aDiv3);
    angle = atan2(sqrt(-q),mbDiv2)*inv3;
    cs = cos(angle);
    sn = sin(angle);
    d[0] = c2Div3 + 2.0*magnitude*cs;
    d[1] = c2Div3 - magnitude*(cs + root3*sn);
    d[2] = c2Div3 - magnitude*(cs - root3*sn);
}   

int fast_eigen3x3(double A[3][3], double V[3][3], double d[3])
{
    const double smallv=1e-12;
    double c0, c1,c2;
    int check;
	double l1, l2, l3;
    double t;
    double a1, a2, a3, b1, b2;
    double da[3];
    check=(absd(A[0][1])<smallv)+(absd(A[0][2])<smallv)+(absd(A[1][2])<smallv);
    if(check>1) { return 0; }

    /* 0 = - det (A - yI) = y^3 + c2 * y^2 + c1 *y + c0 */
    c0 = -(A[0][0]*A[1][1]*A[2][2] + 2*A[0][1]*A[0][2]*A[1][2] - A[0][0] * pow2(A[1][2]) - A[1][1]*pow2(A[0][2]) - A[2][2]*pow2(A[0][1]));
    c1 = A[0][0]*A[1][1] - pow2(A[0][1]) + A[0][0]*A[2][2] -pow2(A[0][2]) + A[1][1]*A[2][2] - pow2(A[1][2]);
    c2 = - (A[0][0] + A[1][1] + A[2][2]);

    /* Solve the roots of  y^3 + c2 * y^2 + c1 *y + c0  */
    roots3(d, c0, c1, c2);

    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);
    /* Sort eigenvalues */
    if(da[0]>=da[1])
    {
        if(da[0]>da[2])
        {
            t=d[0]; d[0]=d[2]; d[2]=t; 
            t=da[0]; da[0]=da[2]; da[2]=t; 
        }
    }
    else if(da[1]>da[2])
    {
        t=d[1]; d[1]=d[2]; d[2]=t; 
        t=da[1]; da[1]=da[2]; da[2]=t; 

    }
    
    if(da[0]>=da[1])
    {
        t=d[0]; d[0]=d[1]; d[1]=t; 
        t=da[0]; da[0]=da[1]; da[1]=t; 
    }

    if((da[1]-da[0])<smallv) { return 0; }
    if((da[2]-da[1])<smallv) { return 0; }
    
    /* Calculate eigen vectors */
    a1=A[0][1]*A[1][2]; a2=A[0][1]*A[0][2]; a3=pow2(A[0][1]);

    b1=A[0][0]-d[0]; b2=A[1][1]-d[0];
    V[0][0]=a1-A[0][2]*b2; V[1][0]=a2-A[1][2]*b1; V[2][0]=b1*b2-a3;

    b1=A[0][0]-d[1]; b2=A[1][1]-d[1];
    V[0][1]=a1-A[0][2]*b2; V[1][1]=a2-A[1][2]*b1; V[2][1]=b1*b2-a3;

    b1=A[0][0]-d[2]; b2=A[1][1]-d[2];
    V[0][2]=a1-A[0][2]*b2; V[1][2]=a2-A[1][2]*b1; V[2][2]=b1*b2-a3;


    /* Eigen vector normalization */
    l1=sqrt(pow2(V[0][0])+ pow2(V[1][0]) + pow2(V[2][0]));
    l2=sqrt(pow2(V[0][1])+ pow2(V[1][1]) + pow2(V[2][1]));
    l3=sqrt(pow2(V[0][2])+ pow2(V[1][2]) + pow2(V[2][2]));

    /* Detect fail : eigenvectors with only zeros */
    if(l1<smallv) { return 0; }
    if(l2<smallv) {	return 0; }
    if(l3<smallv) { return 0; }

    V[0][0]/=l1; V[0][1]/=l2; V[0][2]/=l3;
    V[1][0]/=l1; V[1][1]/=l2; V[1][2]/=l3;
    V[2][0]/=l1; V[2][1]/=l2; V[2][2]/=l3;
    
    /* Succes    */
    return 1;
}

void eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {
    double e[n];
    double da[3];
    double dt, dat;
    double vet[3];
    int i, j;
    
    if(fast_eigen3x3(A, V, d)) { return; }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
    
    /* Sort the eigen values and vectors by abs eigen value */
    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);
    if((da[0]>=da[1])&&(da[0]>da[2]))
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }
    else if((da[1]>=da[0])&&(da[1]>da[2]))  
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
        d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2]; 
    }
    if(da[0]>da[1])
    {
        dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
        d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }
}


void diffusion_scheme_3D_rotation_invariance(float *u,float *u_new, int *dimsu, float *Dxx,float *Dxy,float *Dxz,float *Dyy,float *Dyz,float *Dzz, double dt_d){
	float *j1, *j2, *j3, *ud, *du;
    int npixels;
    float dt=0.5f;
	int i;
    int x, y,z,index;
            
    dt=(float) dt_d;
    npixels=dimsu[0]*dimsu[1]*dimsu[2];
  
 	j1 =mallocf(npixels);
	j2 =mallocf(npixels);
	j3 =mallocf(npixels);
	ud =mallocf(npixels);
	
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
	du =mallocf(npixels);	
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

	for (i=0; i<npixels; i++) 
	{ 
		u_new[i]+=u[i]+du[i]*dt; 
	}

	/* Free memory */
	free(du);
}

#ifdef _WIN32
  unsigned __stdcall   StructureTensor2DiffusionTensorThread(float **Args)  {
#else
  void StructureTensor2DiffusionTensorThread(float **Args)  {
#endif
    /* Matrices of Eigenvector calculation */
    double Ma[3][3];
    double Davec[3][3];
    double Daeig[3];
    
    /* Eigenvector and eigenvalues as scalars */
    double mu1, mu2, mu3, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
    
    /* Amplitudes of diffustion tensor */
    double lambda1, lambda2, lambda3;
      
     /* Denominator  */
    double di;
    
    /* Eps for finite values */
    float eps=(float)1e-20;
    
    /* Loop variable */
    int i;
    
    /* Number of pixels */
    int npixels=1;
    
    float *Jxx, *Jxy, *Jxz, *Jyy, *Jyz, *Jzz;
    float *Dxx, *Dxy, *Dxz, *Dyy, *Dyz, *Dzz;
    int dimsu[3];
    float *dimsu_f, *constants_f, *Nthreads_f, *ThreadID_f;
    double C, m, alpha;
    int ThreadOffset, Nthreads;
    
    Jxx=Args[0];
    Jxy=Args[1];
    Jxz=Args[2];
    Jyy=Args[3];
    Jyz=Args[4];
    Jzz=Args[5];
    Dxx=Args[6];
    Dxy=Args[7];
    Dxz=Args[8];
    Dyy=Args[9];
    Dyz=Args[10];
    Dzz=Args[11];
    dimsu_f=Args[12];
    constants_f=Args[13];
    ThreadID_f=Args[14];
    Nthreads_f=Args[15];
            
    for(i=0;i<3;i++){ dimsu[i]=(int)dimsu_f[i]; }
    C=(double)constants_f[0]; 
    m=(double)constants_f[1]; 
    alpha=(double)constants_f[2];
    ThreadOffset=(int)ThreadID_f[0];
    Nthreads=(int)Nthreads_f[0];    
    
    npixels=dimsu[0]*dimsu[1]*dimsu[2];
    
    for(i=ThreadOffset; i<npixels; i=i+Nthreads) {
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
        di=(mu1-mu3);
        if((di<eps)&&(di>-eps)) { lambda1 = alpha; } else { lambda1 = alpha + (1.0- alpha)*exp(-C/pow(di,(2.0*m))); }
         /*di=(mu2-mu3); */
        /*if((di<eps)&&(di>-eps)) { lambda2 = alpha; } else { lambda2 = alpha + (1.0- alpha)*exp(-C/pow(di,(2.0*m))); } */
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




void StructureTensor2DiffusionTensor(float *Jxx,float *Jxy,float *Jxz,float *Jyy,float *Jyz,float *Jzz,float *Dxx,float *Dxy,float *Dxz,float *Dyy,float *Dyz,float *Dzz,int *dimsu, double C,double m,double alpha) {
    /* ID of Threads */
    float **ThreadID;
    float *ThreadID1;
    float ***ThreadArgs;
    float **ThreadArgs1;
    float Nthreads_f[1]={0};
    float dimsu_f[3];
    float constants_f[3];
    int Nthreads;
    int i;
    
    /* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
	    
    Nthreads=2;
    Nthreads_f[0]=(float)Nthreads;
    for(i=0; i<3; i++) { dimsu_f[i]=(float)dimsu[i]; }
    constants_f[0]=(float)C; constants_f[1]=(float)m; constants_f[2]=(float)alpha;
    
    /* Reserve room for handles of threads in ThreadList  */
	#ifdef _WIN32
		ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
    #else
		ThreadList = (pthread_t*)malloc(Nthreads* sizeof( pthread_t ));
	#endif

    ThreadID = (float **)malloc( Nthreads* sizeof(float *) );
    ThreadArgs = (float ***)malloc( Nthreads* sizeof(float **) );
        
    for (i=0; i<Nthreads; i++) {
        /*  Make Thread ID  */
        ThreadID1= (float *)malloc( 1* sizeof(float) );
        ThreadID1[0]=(float)i;
        ThreadID[i]=ThreadID1;
        
        /*  Make Thread Structure  */
        ThreadArgs1 = (float **)malloc( 16* sizeof( float * ) );
        ThreadArgs1[0]=Jxx;
        ThreadArgs1[1]=Jxy;
        ThreadArgs1[2]=Jxz;
        ThreadArgs1[3]=Jyy;
        ThreadArgs1[4]=Jyz;
        ThreadArgs1[5]=Jzz;
        ThreadArgs1[6]=Dxx;
        ThreadArgs1[7]=Dxy;
        ThreadArgs1[8]=Dxz;
        ThreadArgs1[9]=Dyy;
        ThreadArgs1[10]=Dyz;
        ThreadArgs1[11]=Dzz;
        ThreadArgs1[12]=dimsu_f;
        ThreadArgs1[13]=constants_f;
        ThreadArgs1[14]=ThreadID[i];
        ThreadArgs1[15]=Nthreads_f;
       
        /* Start a Thread  */
        ThreadArgs[i]=ThreadArgs1;
		#ifdef _WIN32
			ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &StructureTensor2DiffusionTensorThread, ThreadArgs[i] , 0, NULL );
		#else
			pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &StructureTensor2DiffusionTensorThread, ThreadArgs[i]);
		#endif
    }
    
    #ifdef _WIN32
		for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
		for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
	#else
		for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
	#endif
    
    for (i=0; i<Nthreads; i++) {
        free(ThreadArgs[i]);
        free(ThreadID[i]);
    }
    
    free(ThreadArgs);
    free(ThreadID );
    free(ThreadList);

    
}


