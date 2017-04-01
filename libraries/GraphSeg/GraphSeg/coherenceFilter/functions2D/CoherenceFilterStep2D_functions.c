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

        
__inline double pow2(double a) { return a*a;}

void imfilter1D_double(double *I, int lengthI, double *H, int lengthH, double *J) {
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

void imfilter2D_double(double *I, int * sizeI, double *H, int lengthH, double *J) {
    int y, x, i, y2;
    double *Irow, *Crow;
    int index=0, line=0;
    double *RCache;
    int *nCache;
    int hks, offset, offset2;
    RCache=(double *)malloc(lengthH*sizeI[0]*sizeof(double));
    for(i=0; i<lengthH*sizeI[0]; i++) { RCache[i]=0; }
    nCache=(int *)malloc(lengthH*sizeof(int));
    for(i=0; i<lengthH; i++) { nCache[i]=0; }
    hks=((lengthH-1)/2);
    for(y=0; y<min(hks,sizeI[1]); y++) {
        Irow=&I[index];
        Crow=&RCache[line*sizeI[0]];
        imfilter1D_double(Irow, sizeI[0], H, lengthH, Crow);
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
        imfilter1D_double(Irow, sizeI[0], H, lengthH, Crow);
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
        imfilter1D_double(Irow, sizeI[0], H, lengthH, Crow);
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

void imfilter2Dcolor_double(double *I, int * sizeI, double *H, int lengthH, double *J) {
	int i, index;
	for(i=0; i<sizeI[2]; i++)
	{
		index=i*(sizeI[0]*sizeI[1]);
		imfilter2D_double(&I[index], sizeI, H, lengthH, &J[index]);
	}
}

void GaussianFiltering2Dcolor_double(double *I, double *J, int *dimsI, double sigma, double kernel_size)
{
	int kernel_length,i;
    double x, *H, totalH=0;
	
	/* Construct the 1D gaussian kernel */
	if(kernel_size<1) { kernel_size=1; }
    kernel_length=(int)(2*ceil(kernel_size/2)+1);
	H = (double *)malloc(kernel_length*sizeof(double));
	x=-ceil(kernel_size/2);
	for (i=0; i<kernel_length; i++) { H[i]=exp(-((x*x)/(2*(sigma*sigma)))); totalH+=H[i]; x++; }
	for (i=0; i<kernel_length; i++) { H[i]/=totalH; }
	
	/* Do the filtering */
	imfilter2Dcolor_double(I, dimsI, H, kernel_length, J);
    /* Clear memory gaussian kernel */
	free(H);
}

void GaussianFiltering2D_double(double *I, double *J, int *dimsI, double sigma, double kernel_size)
{
	int kernel_length,i;
    double x, *H, totalH=0;
	
	/* Construct the 1D gaussian kernel */
	if(kernel_size<1) { kernel_size=1; }
    kernel_length=(int)(2*ceil(kernel_size/2)+1);
	H = (double *)malloc(kernel_length*sizeof(double));
	x=-ceil(kernel_size/2);
	for (i=0; i<kernel_length; i++) { H[i]=exp(-((x*x)/(2*(sigma*sigma)))); totalH+=H[i]; x++; }
	for (i=0; i<kernel_length; i++) { H[i]/=totalH; }
	
	/* Do the filtering */
	imfilter2D_double(I, dimsI, H, kernel_length, J);
    /* Clear memory gaussian kernel */
	free(H);
}

double * mallocd(int a)
{
    double *ptr = malloc(a*sizeof(double));
    if (ptr == NULL) { mexErrMsgTxt("Out of Memory"); }
    return ptr;
}



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

void StructureTensor2D(double *ux, double  *uy,  double **J, int *dimsu, double rho) {
    int npixels2;
    int i,j,offset;
    double *Jxx, *Jxy, *Jyy; /* Structure tensor */
       
    npixels2=dimsu[0]*dimsu[1];
    
    Jxx = mallocd(npixels2);  
    Jyy = mallocd(npixels2); 
    Jxy = mallocd(npixels2);  
    
    /* J(grad u_sigma) */
	for(i=0; i<npixels2; i++)
	{
		Jxx[i] = ux[i]*ux[i]; Jxy[i] = ux[i]*uy[i]; Jyy[i] = uy[i]*uy[i];
	}
	for(j=1; j<dimsu[2]; j++)
	{
		offset=j*npixels2;
		for(i=0; i<npixels2; i++)
		{
			Jxx[i] += ux[i+offset]*ux[i+offset]; 
			Jxy[i] += ux[i+offset]*uy[i+offset];
			Jyy[i] += uy[i+offset]*uy[i+offset];
		}
	}
	
    /* Do the gaussian smoothing */
    J[0]=mallocd(npixels2); 
    GaussianFiltering2D_double(Jxx, J[0], dimsu, rho, 4*rho);
    free(Jxx);
    J[1]=mallocd(npixels2); 
    GaussianFiltering2D_double(Jyy, J[1], dimsu, rho, 4*rho);
    free(Jyy);
    J[2]=mallocd(npixels2); 
    GaussianFiltering2D_double(Jxy, J[2], dimsu, rho, 4*rho);
    free(Jxy);
 }
 
 
void EigenVector2D(double Jxx, double Jxy, double Jyy, double *Result)
{
	double mu1, mu2, v2x, v2y, v1x, v1y;
	double tmp, mag;
    
	/* Compute the eigenvectors of J, v1 and v2 */
	tmp = sqrt(pow2(Jxx - Jyy) + 4*pow2(Jxy));
	v2x = 2*Jxy; v2y = Jyy - Jxx + tmp;

	/* Normalize */
	mag = sqrt(pow2(v2x) + pow2(v2y));
	if(mag!=0) { v2x = v2x/mag; v2y = v2y/mag; }
	
	/* The eigenvectors are orthogonal */
	v1x = -v2y; v1y = v2x;

	/* Compute the eigenvalues */
	mu1 = 0.5*(Jxx + Jyy + tmp);
	mu2 = 0.5*(Jxx + Jyy - tmp);
	
	/* Make output structure */
	Result[0]=mu2; Result[1]=mu1; Result[2]=v2x; Result[3]=v2y; Result[4]=v1x; Result[5]=v1y;
}
	 
void StructureTensor2DiffusionTensor(double *Jxx,double *Jxy,double *Jyy,double *Dxx,double *Dxy,double *Dyy,int *dimsu,double C, double m, double alpha)
{
	double mu1, mu2, v2x, v2y, v1x, v1y;
	double Result[6];
    int i;
    double di;
    double lambda1, lambda2;
    /* Eps for finite values */
    double eps=(double)1e-20;
    
	for (i=0; i<(dimsu[0]*dimsu[1]); i++) 
	{
		EigenVector2D(Jxx[i], Jxy[i], Jyy[i], Result);
		mu2=Result[0]; mu1=Result[1]; v2x=Result[2]; v2y=Result[3]; v1x=Result[4]; v1y=Result[5];
		
		/* Scaling of diffusion tensor */
        di=(mu1-mu2);
        if((di<eps)&&(di>-eps)) { lambda1 = alpha; } else { lambda1 = alpha + (1.0- alpha)*exp(-C/pow(di,(2.0*m))); }
        lambda2 = alpha;

        /* Construct the diffusion tensor */
        Dxx[i] = lambda1*v1x*v1x + lambda2*v2x*v2x;
        Dxy[i] = lambda1*v1x*v1y + lambda2*v2x*v2y;
		Dyy[i] = lambda1*v1y*v1y + lambda2*v2y*v2y;
        
	}
}


void diffusion_scheme_2D_rotation_invariance(double *u,double *u_new,int *dimsu,double *Dxx,double *Dxy,double *Dyy, double dt)
{
	double  *j1, *j2, *ud, *du;
    int npixels2, npixels3;
    int i,j, offset, offset2;
    int x, y;
            
    npixels2=dimsu[0]*dimsu[1];
    npixels3=dimsu[0]*dimsu[1]*dimsu[2];
  
 	j1 =mallocd(npixels3);
	j2 =mallocd(npixels3);
	ud =mallocd(npixels3);
	
	/* 3 : Calculate the flux components */
	/* j1 = Dxx .* ux + Dxy .*uy; */
	/* j2 = Dxy .* ux + Dyy .*uy; */
	
	gradient2Dx_color(u, dimsu, ud);
	for(j=0; j<dimsu[2]; j++)
	{
		offset=j*npixels2;
		for (i=0; i<npixels2; i++) 
		{ 
			j1[i+offset]=Dxx[i]*ud[i+offset]; j2[i+offset]=Dxy[i]*ud[i+offset]; 
		}
	}
	
	gradient2Dy_color(u, dimsu, ud);
	for(j=0; j<dimsu[2]; j++)
	{
		offset=j*npixels2;
		for (i=0; i<npixels2; i++) 
		{ 
			j1[i+offset]+=Dxy[i]*ud[i+offset]; j2[i+offset]+=Dyy[i]*ud[i+offset]; 
		}
	}
    
	/* j1(1,:)=0; j1(end,:)=0; j1(:,1)=0; j1(:,end)=0; */
	/* j2(1,:)=0; j2(end,:)=0; j2(:,1)=0; j2(:,end)=0; */
    
	for(j=0; j<dimsu[2]; j++)
	{
		offset2=j*npixels2;
		for (x=0; x<dimsu[0]; x++) 
		{ 
			j1[x+offset2]=0; j2[x+offset2]=0;
		}
		offset=(dimsu[1]-1)*dimsu[0];
		for (x=0; x<dimsu[0]; x++) 
		{ 
			j1[x+offset+offset2]=0; j2[x+offset+offset2]=0;
		}
		offset=0; 
		for (y=0; y<dimsu[1]; y++) 
		{ 
			j1[offset+offset2]=0; j2[offset+offset2]=0;
			j1[offset+dimsu[0]-1+offset2]=0; j2[offset+dimsu[0]-1+offset2]=0;
			offset+=dimsu[0];
		}
    }
    
    /* 4 : Calculate ... by means of the optimized derivative filters */
	/* du = derivatives(j1,'x')+derivatives(j2,'y'); */
	du =mallocd(npixels3);	
	gradient2Dx_color(j1, dimsu, du);
	gradient2Dy_color(j2, dimsu, ud);
	for (i=0; i<npixels3; i++) { du[i]+=ud[i]; }
	
    /* Free memory */
	free(j1);
	free(j2);
	free(ud);
	
	/* 5 : Update in an explicit way. */
	/* u=u+du*dt; */

	for (i=0; i<npixels3; i++) 
	{ 
		u_new[i]+=u[i]+du[i]*dt; 
	}

	/* Free memory */
	free(du);
}
 
