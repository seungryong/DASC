#include <mex.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#define max(X,Y) ((X) > (Y) ? (X) : (Y)) 
#define min(X,Y) ((X) < (Y) ? (X) : (Y))  

int num_iterations;

static inline void domaintransform_runfilter(float *img, float *V_dHdx, float *V_dVdy, float *img_out, int height, int width);
static inline void TransformedDomainRecursiveFilter_Horizontal(float *I, float *V, int iter, int height, int width);
static inline void diff(float *img, float *img_out, int dim, int height, int width);
static inline void image_transpose(float *img, float *img_out, int height, int width);

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *image;
    double *fVol;
    double *rp1, *rp2;
    float dem_left, dem_right;
    float num_corrSurf, dem_corrSurf;
    float fnorm;
    float sigma_s, sigma_r;
    int i, j, s;
    int x, y;
    int f_dim;
    int m, n, ind, ind1, i1, j1;
    int height, width;
    
    image = mxGetPr(prhs[0]); 
	rp1 = mxGetPr(prhs[1]);
    rp2 = mxGetPr(prhs[2]);  
    sigma_s = (float)(*mxGetPr(prhs[3]));
    sigma_r = (float)(*mxGetPr(prhs[4]));
    num_iterations = (int)(*mxGetPr(prhs[5]));
    
	const mwSize *Size = mxGetDimensions(prhs[0]);
    height = Size[0];
    width = Size[1];
    f_dim = (int)mxGetM(prhs[1]); 
    
    mwSize dimK = 3;      
	const mwSize dims_fVol[] = {height,width,f_dim};
    plhs[0] = mxCreateNumericArray(dimK,dims_fVol,mxDOUBLE_CLASS,mxREAL);    
    fVol = mxGetPr(plhs[0]);  
    
    int size_image = (int)height*width*sizeof(float);
    
    float *I = (float*)malloc(size_image);
    float *II = (float*)malloc(size_image);
    float *I_adaptive_mean = (float*)malloc(size_image);
    float *II_adaptive_mean = (float*)malloc(size_image);
    float *J_adaptive_mean = (float*)malloc(size_image);
    float *JJ_adaptive_mean = (float*)malloc(size_image);
    float *IJ_adaptive_mean = (float*)malloc(size_image);
    float *J = (float*)malloc(size_image);
    float *JJ = (float*)malloc(size_image);
    float *IJ = (float*)malloc(size_image);
    float *fout = (float*)malloc(height*width*f_dim*sizeof(float));    
        
    int *diff_rp = (int*)malloc(f_dim*2*sizeof(int));
    
    memset(I,0,size_image);
    memset(II,0,size_image);
    memset(I_adaptive_mean,0,size_image);
    memset(II_adaptive_mean,0,size_image);
    memset(J_adaptive_mean,0,size_image);
    memset(JJ_adaptive_mean,0,size_image);
    memset(IJ_adaptive_mean,0,size_image);
    memset(J,0,size_image);
    memset(JJ,0,size_image);
    memset(IJ,0,size_image);
    memset(fout,0,height*width*f_dim*sizeof(float));
    
    float *dIcdx = (float*)malloc(size_image);
    float *dIcdy = (float*)malloc(size_image);
    memset(dIcdx,0,size_image);
    memset(dIcdy,0,size_image);
    
    float *dHdx = (float*)malloc(size_image);
    float *dVdy = (float*)malloc(size_image);
    memset(dHdx,0,size_image);
    memset(dVdy,0,size_image);

    for( i=0;i<height*width;i++ )
    {
        I[i] = image[i];
        II[i] = I[i]*I[i];
    }
    
	for (s = 0; s<f_dim; s++)
	{
		diff_rp[s] = rp2[s] - rp1[s];
		diff_rp[s + f_dim] = rp2[s + f_dim] - rp1[s + f_dim];
	}
    
    // domain transform pre-computation
	diff(I,dIcdx,2,height,width);    
    diff(I,dIcdy,1,height,width);    
    
	for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
        {
            dHdx[y + height*x] = 1.f + sigma_s/sigma_r * abs(dIcdx[y + height*x]);
            dVdy[x + width*y] = 1.f + sigma_s/sigma_r * abs(dIcdy[y + height*x]);
        }           
    
    float *V_dHdx = (float*)malloc((int)height*width*num_iterations*sizeof(float));
    float *V_dVdy = (float*)malloc((int)height*width*num_iterations*sizeof(float));
    memset(V_dHdx,0,(int)height*width*num_iterations*sizeof(float));
    memset(V_dVdy,0,(int)height*width*num_iterations*sizeof(float));
    
    int N = num_iterations;    
    float sigma_H = sigma_s;
        
    for (i = 0; i < num_iterations; i++)
    {
        float sigma_H_i = sigma_H * sqrt(3.f) * (float)pow( 2.f, (float)(N - (i+1)) ) / sqrt( pow(4.f, (float)N ) - 1 );
        float a = exp( - sqrt(2.f) / sigma_H_i );
        
        for (y = 0; y < height; y++)
            for (x = 0; x < width; x++)
            {
                V_dHdx[y + height*(x + width*i)] = (float)pow(a,(float)dHdx[y + height*x]);
                V_dVdy[x + width*(y + height*i)] = (float)pow(a,(float)dVdy[x + width*y]);
            }        
    }  

    domaintransform_runfilter(I,V_dHdx,V_dVdy,I_adaptive_mean,height,width);
    domaintransform_runfilter(II,V_dHdx,V_dVdy,II_adaptive_mean,height,width);
    
	for (s = 0; s<f_dim; s++)
	{
		m = diff_rp[s];
		n = diff_rp[s + f_dim];

		for (i = 0; i<height; i++)
            for (j = 0; j<width; j++)
                if (i + m > -1 && i + m < height && j + n > -1 && j + n < width)
                {
                    ind = i + height*j;
                    ind1 = (i + m) + height*(j + n);
                    J[ind] = I[ind1];
                    JJ[ind] = I[ind1]*I[ind1];
                    IJ[ind] = I[ind]*I[ind1];
                }

        domaintransform_runfilter(J,V_dHdx,V_dVdy,J_adaptive_mean,height,width);
        domaintransform_runfilter(JJ,V_dHdx,V_dVdy,JJ_adaptive_mean,height,width);
        domaintransform_runfilter(IJ,V_dHdx,V_dVdy,IJ_adaptive_mean,height,width);

		for (i = 0; i<height; i++)
		for (j = 0; j<width; j++)
		{
			i1 = (int)(i + rp1[s]);
			j1 = (int)(j + rp1[s + f_dim]);
            
			if (i1 > 0 && i1 < height && j1 > 0 && j1 < width)
            {
				ind = i1 + height*j1;
				num_corrSurf = IJ_adaptive_mean[ind] - I_adaptive_mean[ind] * J_adaptive_mean[ind];
				dem_left = II_adaptive_mean[ind] - I_adaptive_mean[ind] * I_adaptive_mean[ind];
				dem_right = JJ_adaptive_mean[ind] - J_adaptive_mean[ind] * J_adaptive_mean[ind];
                dem_corrSurf = sqrt(dem_left*dem_right);
                if( dem_corrSurf>1e-10 )
                    fout[i + height*(j + width*s)] = exp(-(1-num_corrSurf/dem_corrSurf)/0.5);
                else
                    fout[i + height*(j + width*s)] = 1;                
                
			}
		}
	}

    float sqrt_fnorm;
	for (i = 0; i<height; i++)
        for (j = 0; j<width; j++)
        {
            fnorm = 0;
            for (s = 0; s<f_dim; s++)
                fnorm += fout[i + height*(j + width*s)] * fout[i + height*(j + width*s)];
            sqrt_fnorm = sqrt(fnorm);
            
            if( sqrt_fnorm > 1e-10 )
                for (s = 0; s<f_dim; s++)
                    fVol[i + height*(j + width*s)] = fout[i + height*(j + width*s)] / sqrt_fnorm;
        }

	free(I);
	free(II);
	free(I_adaptive_mean);
	free(II_adaptive_mean);
	free(J_adaptive_mean);
	free(JJ_adaptive_mean);
	free(IJ_adaptive_mean);
	free(J);
	free(JJ);
	free(IJ);
    free(fout);
    
    free(dIcdx); 
    free(dIcdy); 
    free(dHdx); 
    free(dVdy); 
    free(V_dHdx); 
    free(V_dVdy); 
    
	free(diff_rp);
}

static inline void domaintransform_runfilter(float *img, float *V_dHdx, float *V_dVdy, float *img_out, int height, int width)
{
    int size_image = (int)height*width*sizeof(float);
    float *F = (float*)malloc(size_image);
    float *FT = (float*)malloc(size_image);
    memset(F,0,size_image);
	memset(FT,0,size_image);
    memset(img_out,0,size_image);
    
    for (int i = 0; i < height*width; i++)
        F[i] = img[i];  
    
    // Perform the filtering    
    for (int iter = 0; iter < num_iterations; iter++)
    {
        TransformedDomainRecursiveFilter_Horizontal(F,V_dHdx,iter,height,width);
        image_transpose(F,FT,height,width);
        TransformedDomainRecursiveFilter_Horizontal(FT,V_dVdy,iter,width,height);
        image_transpose(FT,F,width,height);
    }      
    
    for (int i = 0; i < height*width; i++)
        img_out[i] = F[i];  
    
    free(F);
    free(FT);
}

static inline void TransformedDomainRecursiveFilter_Horizontal(float *I, float *V, int iter, int height, int width)
{
    // Left -> Right filter.
	for(int y=0;y<height;y++)
        for(int x=1;x<width;x++)
            I[y + height*x] = I[y + height*x] + V[y + height*(x+width*iter)] * ( I[y + height*(x-1)] - I[y + height*x] );
    
    // Right -> Left filter.
	for(int y=0;y<height;y++)
        for(int x=width-2;x>=0;x--)
            I[y + height*x] = I[y + height*x] + V[y + height*((x+1)+width*iter)] * ( I[y + height*(x+1)] - I[y + height*x] );   
}

static inline void diff(float *img, float *img_out,  int dim, int height, int width)
{
    int dx, dy;
    
    if (dim==1){
        dy = 1;
        dx = 0;
    }
    else{
        dy = 0;
        dx = 1;
    }    
    
    for (int y = dy; y < height; y++)
        for (int x = dx; x < width; x++)
            img_out[y + height*x] = img[(y-dy) + height*(x-dx)] - img[y + height*x];   
}

static inline void image_transpose(float *img, float *img_out, int height, int width)
{
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            img_out[x + width*y] = img[y + height*x];   
}