#include "mex.h"
#include "math.h"
#ifndef min
#define min(a,b)        ((a) < (b) ? (a): (b))
#endif
#ifndef max
#define max(a,b)        ((a) > (b) ? (a): (b))
#endif
#define absd(a) ((a)>(-a)?(a):(-a))
#define mwSize int

/* Stencil= [3 10 3; 0 0 0; -3 -10 -3]/32;   */
/* The matlab mex function */



void gradient2Dx(double *I, int *sizeI, double *Ix)
{
    int x,y;
    int yp, yn;
    int indexn, indexc, indexp;
    double *Irow;
    const double smoothfilter[2]={0.093750,0.31250};
    /*const double smoothfilter[2]={0.07322330470336310700,0.35355339059327379000}; */
    Irow=(double *)malloc(sizeI[0]*sizeof(double));
    
    for(y=0; y<sizeI[1]; y++)
    {
      /* Smooth in y - direction */
        yn=max(y-1,0);
        yp=min(y+1,sizeI[1]-1);
        indexn=yn*sizeI[0]; indexc=y*sizeI[0]; indexp=yp*sizeI[0];
        for(x=0; x<sizeI[0]; x++)
        {
            Irow[x] =smoothfilter[0]*I[indexn+x];
            Irow[x]+=smoothfilter[1]*I[indexc+x];
            Irow[x]+=smoothfilter[0]*I[indexp+x];
        }
      /* Gradient in x - direction */
        Ix[indexc]=2*(Irow[1]-Irow[0]);
        for(x=1; x<(sizeI[0]-1); x++)
        {
            Ix[indexc+x]=Irow[x+1]-Irow[x-1];
        }
        Ix[indexc+sizeI[0]-1]=2*(Irow[sizeI[0]-1]-Irow[sizeI[0]-2]);
    }
    free(Irow);
}

void gradient2Dy(double *I, int *sizeI, double *Iy)
{
    int x,y;
    int indexy, indexyy, indexn, indexc, indexp;
    double *Irow;
    int row_select=0, row_select_p1=0, row_select_p2=0;
    const double smoothfilter[2]={0.093750,0.31250};
    /*const double smoothfilter[2]={0.07322330470336310700,0.35355339059327379000}; */
    Irow=(double *)malloc(sizeI[0]*4*sizeof(double));
    
    for(y=0; y<sizeI[1]; y++)
    {
        /* Smooth in x - direction */
        indexyy=y*sizeI[0];
        indexy=row_select*sizeI[0];
        Irow[indexy]=I[indexyy]*(smoothfilter[0]+smoothfilter[1])+I[indexyy+1]*smoothfilter[0];
        for(x=1; x<(sizeI[0]-1); x++)
        {
            indexy++; indexyy++;
            Irow[indexy] =smoothfilter[0]*I[indexyy-1];
            Irow[indexy]+=smoothfilter[1]*I[indexyy];
            Irow[indexy]+=smoothfilter[0]*I[indexyy+1];
        }
        indexy++; indexyy++;
        Irow[indexy]=I[indexyy]*(smoothfilter[0]+smoothfilter[1])+I[indexyy-1]*smoothfilter[0];
       
        /* Gradient in y - direction */
        if(y==1) /* Forward */        
        {
            for(x=0; x<sizeI[0]; x++) { Iy[x]=2*(Irow[x+sizeI[0]]-Irow[x]); }
        }
        else if(y>1) /* Central */
        {
            indexn=row_select_p2*sizeI[0]; indexp=row_select*sizeI[0]; indexyy=(y-1)*sizeI[0];
            for(x=0; x<sizeI[0]; x++) { Iy[x+indexyy]=Irow[x+indexp]-Irow[x+indexn];  }
        }
        if(y==(sizeI[1]-1)) /* Backward */
        {
            indexn=row_select_p1*sizeI[0]; indexc=row_select*sizeI[0]; indexyy=(sizeI[1]-1)*sizeI[0];
            for(x=0; x<sizeI[0]; x++)  { Iy[x+indexyy]=2*(Irow[x+indexc]-Irow[x+indexn]);  }
        }
        row_select_p2=row_select_p1; row_select_p1=row_select; row_select++; if(row_select>3) { row_select=0; }
    }
    free(Irow);
}

void gradient3Dx_double(double *I, int *sizeI, double *Ix)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    double *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const double smoothfilter[3]={0.187500,0.625000,0.187500};
    const double derivafilter[3]={-0.5,0,0.5};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(double *)malloc(4*nSlice*sizeof(double));
    Irow=(double *)malloc(sizeI[0]*sizeof(double));
    
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

void gradient3Dy_double(double *I, int *sizeI, double *Iy)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    double *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const double smoothfilter[3]={0.187500,0.625000,0.187500};
    const double derivafilter[3]={-0.5,0,0.5};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(double *)malloc(4*nSlice*sizeof(double));
    Irow=(double *)malloc(sizeI[0]*sizeof(double));
    
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

void gradient3Dz_double(double *I, int *sizeI, double *Iz)
{
    int x,y,z;
    int xp, xn, yp, yn;
    int i;
    int indexn, indexc, indexp;
    double *Irow, *Islices;
    int nSlice;
    int offsetz,offset_slice;
    int slice_select=0, slice_select_p1=0, slice_select_p2=0;
    
    const double smoothfilter[3]={0.187500,0.625000,0.187500};
    const double derivafilter[3]={-0.5,0,0.5};
    
    nSlice=sizeI[0]*sizeI[1];
    Islices=(double *)malloc(4*nSlice*sizeof(double));
    Irow=(double *)malloc(sizeI[0]*sizeof(double));
    
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

void gradient2Dx_color(double *I, int *sizeI, double *J)
{
    int i, index;
	for(i=0; i<sizeI[2]; i++)
	{
		index=i*(sizeI[0]*sizeI[1]);
		gradient2Dx(&I[index], sizeI, &J[index]);
    }
}

void gradient2Dy_color(double *I, int *sizeI, double *J)
{
    int i, index;
	for(i=0; i<sizeI[2]; i++)
	{
		index=i*(sizeI[0]*sizeI[1]);
		gradient2Dy(&I[index], sizeI, &J[index]);
    }
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{
    double *I_double, *J_double;
    float *I_float, *J_float;
    char *option;
    const mwSize *dimsI_const;
    int dimsI[3], ndimsI;

    /* Get image input size */
    ndimsI=mxGetNumberOfDimensions(prhs[0]);
    if((ndimsI<2)||(ndimsI>3)) { mexErrMsgTxt("Image must be 2D"); }
    dimsI_const = mxGetDimensions(prhs[0]);
	dimsI[0]=dimsI_const[0];
    dimsI[1]=dimsI_const[1];
    if(ndimsI>2) { dimsI[2]=dimsI_const[2]; } else { dimsI[2]=1; }
    
    /* Assign pointers to input. */
    option=(char*)mxGetData(prhs[1]);
        
    
    if(mxGetClassID(prhs[0])==mxDOUBLE_CLASS) 
    {       
        /* Create output matrix */
        plhs[0] = mxCreateNumericArray(ndimsI, dimsI, mxDOUBLE_CLASS, mxREAL);
        I_double=(double*)mxGetData(prhs[0]);
        J_double=(double*)mxGetData(plhs[0]);
 
        switch (option[0])
        {
         case 'x': case 'X':
             if(dimsI[2]>3) { gradient3Dx_double(I_double, dimsI, J_double); } else { gradient2Dx_color(I_double, dimsI, J_double); }
           break;
         case 'y': case 'Y':
             if(dimsI[2]>3) { gradient3Dy_double(I_double, dimsI, J_double); } else { gradient2Dy_color(I_double, dimsI, J_double); }
           break;
         case 'z': case 'Z':
             if(dimsI[2]>3) { gradient3Dz_double(I_double, dimsI, J_double); } 
           break;
         default:
           printf("\n  option not defined");
           break;
        }
    }
    else if(mxGetClassID(prhs[0])==mxSINGLE_CLASS) 
    {
        /* Create output matrix */
        plhs[0] = mxCreateNumericArray(ndimsI, dimsI, mxSINGLE_CLASS, mxREAL);
        I_float=(float*)mxGetData(prhs[0]);
        J_float=(float*)mxGetData(plhs[0] );
        
        switch (option[0])
        {
         case 'x': case 'X':
             if(dimsI[2]>3) { gradient3Dx_float(I_float, dimsI, J_float); } 
           break;
         case 'y': case 'Y':
             if(dimsI[2]>3) { gradient3Dy_float(I_float, dimsI, J_float); } 
           break;
         case 'z': case 'Z':
             if(dimsI[2]>3) { gradient3Dz_float(I_float, dimsI, J_float); } 
           break;
         default:
           printf("\n  option not defined");
           break;
        }
    }
     
}
        
