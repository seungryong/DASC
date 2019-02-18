#include <mex.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#define max(X,Y) ((X) > (Y) ? (X) : (Y)) 
#define min(X,Y) ((X) < (Y) ? (X) : (Y))  

typedef struct Param{
	int row, col, r, len, downsize;
    double epsil;
} PARAM;

typedef struct Precompute{
    double *integral_img;
    double *glob_I_sub;
    double *glob_I_mean;
    double *glob_II_mean;
    double *glob_I_var;
} PRECOMPUTE;

typedef struct Runfilter{
    double *a_mean; 
    double *b_mean; 
    double *p_mean; 
    double *Ip_mean;     
    double *p_sub; 
    double *a_mean_large; 
    double *b_mean_large; 
} RUNFILTER;

static inline void imresize_nearest(double *img, double *img_out, PARAM *param1, PARAM *param2);
static inline void imresize_bilinear(double *img, double *img_out, PARAM *param1, PARAM *param2);
static inline void boxfilter(double *img, double *img_out, double *integral_img, PARAM *param_sub);
static inline void guidedfilter_precompute(double *I, PRECOMPUTE *preco, PARAM *param, PARAM *param_sub);
static inline void guidedfilter_runfilter(double *I, double *p, double *I_out, PRECOMPUTE *preco, RUNFILTER *runf, PARAM *param, PARAM *param_sub);

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *image;
    double *fVol;
    double *rp1, *rp2;
    double dem_left, dem_right;
    double num_corrSurf, dem_corrSurf;
    double fnorm;
    int max_radius;
    int row_o, col_o;
    int i, j, m, n, s, f_dim, i1, j1, i2, j2, ii, jj;
    int ind, ind1;
    
    PARAM *param = (PARAM*)malloc(sizeof(PARAM));
  
    image = mxGetPr(prhs[0]);
    max_radius = (int)(*mxGetPr(prhs[1]));
    param->r = (int)(*mxGetPr(prhs[2]));
	rp1 = mxGetPr(prhs[3]);
    rp2 = mxGetPr(prhs[4]);  
    param->epsil = (double)(*mxGetPr(prhs[5]));
    param->downsize = (int)(*mxGetPr(prhs[6]));
    
	const mwSize *Size = mxGetDimensions(prhs[0]);
    row_o = Size[0];
    col_o = Size[1];
    
    if(row_o%2 == 1) param->row = row_o-1;
    else param->row = row_o;
    if(col_o%2 == 1) param->col = col_o-1;
    else param->col = col_o;
    
    param->len = param->row*param->col;
    f_dim = (int)mxGetM(prhs[3]); 
    
    mwSize dimK = 3;      
	const mwSize dims_fVol[] = {row_o,col_o,f_dim};
    plhs[0] = mxCreateNumericArray(dimK,dims_fVol,mxDOUBLE_CLASS,mxREAL);    
    fVol = mxGetPr(plhs[0]);  
    
    printf("DASC description...\n");
   
	PARAM *param_sub = (PARAM*)malloc(sizeof(PARAM));
    param_sub->row = (int)param->row/param->downsize;
    param_sub->col = (int)param->col/param->downsize;   
	param_sub->r = (int)param->r/param->downsize;
    param_sub->len = param_sub->row*param_sub->col;
    
    int size_param = (int)param->row*param->col*sizeof(double);
    int size_param_sub = (int)param_sub->row*param_sub->col*sizeof(double);
    
    double *I = (double*)malloc(size_param);
    double *II = (double*)malloc(size_param);
    double *I_adaptive_mean = (double*)malloc(size_param);
    double *II_adaptive_mean = (double*)malloc(size_param);
    double *J_adaptive_mean = (double*)malloc(size_param);
    double *JJ_adaptive_mean = (double*)malloc(size_param);
    double *IJ_adaptive_mean = (double*)malloc(size_param);
    double *J = (double*)malloc(size_param);
    double *JJ = (double*)malloc(size_param);
    double *IJ = (double*)malloc(size_param);
    double *fout = (double*)malloc(param->row*param->col*f_dim*sizeof(double));
    
    PRECOMPUTE *preco = (PRECOMPUTE*)malloc(sizeof(PRECOMPUTE));    
    preco->integral_img = (double*)malloc(size_param_sub);
    preco->glob_I_sub = (double*)malloc(size_param_sub);    
	preco->glob_I_mean = (double*)malloc(size_param_sub);
	preco->glob_II_mean = (double*)malloc(size_param_sub);
	preco->glob_I_var = (double*)malloc(size_param_sub);
    
    RUNFILTER *runf = (RUNFILTER*)malloc(sizeof(RUNFILTER));   
	runf->a_mean = (double*)malloc(size_param_sub); 
    runf->b_mean = (double*)malloc(size_param_sub); 
    runf->p_mean = (double*)malloc(size_param_sub); 
    runf->Ip_mean = (double*)malloc(size_param_sub);     
    runf->p_sub = (double*)malloc(size_param_sub);       
    runf->a_mean_large = (double*)malloc(size_param); 
    runf->b_mean_large = (double*)malloc(size_param); 
    
    int *diff_rp = (int*)malloc(f_dim*2*sizeof(int));
    
    memset(I,0,size_param);
    memset(II,0,size_param);
    memset(I_adaptive_mean,0,size_param);
    memset(II_adaptive_mean,0,size_param);
    memset(J_adaptive_mean,0,size_param);
    memset(JJ_adaptive_mean,0,size_param);
    memset(IJ_adaptive_mean,0,size_param);
    memset(J,0,size_param);
    memset(JJ,0,size_param);
    memset(IJ,0,size_param);
    memset(fout,0,param->row*param->col*f_dim*sizeof(double));

	for (i = 0; i<param->row; i++)
	for (j = 0; j<param->col; j++)
        I[i+param->row*j] = image[i+row_o*j];

    for( i=0;i<param->len;i++ )
        II[i] = I[i]*I[i];
    
	for (s = 0; s<f_dim; s++)
	{
		diff_rp[s] = rp2[s] - rp1[s];
		diff_rp[s + f_dim] = rp2[s + f_dim] - rp1[s + f_dim];
	}
       
	guidedfilter_precompute(I, preco, param, param_sub);
	guidedfilter_runfilter(I, I, I_adaptive_mean, preco, runf, param, param_sub);
	guidedfilter_runfilter(I, II, II_adaptive_mean, preco, runf, param, param_sub);
    
	for (s = 0; s<f_dim; s++)
	{
		m = diff_rp[s];
		n = diff_rp[s + f_dim];

		for (i = 0; i<param->row; i++)
		for (j = 0; j<param->col; j++)
		if (i + m > -1 && i + m < param->row && j + n > -1 && j + n < param->col)
		{
			ind = i + param->row*j;
			ind1 = (i + m) + param->row*(j + n);
			J[ind] = I[ind1];
			JJ[ind] = I[ind1]*I[ind1];
			IJ[ind] = I[ind]*I[ind1];
		}

		guidedfilter_runfilter(I, J, J_adaptive_mean, preco, runf, param, param_sub);
		guidedfilter_runfilter(I, JJ, JJ_adaptive_mean, preco, runf, param, param_sub);
		guidedfilter_runfilter(I, IJ, IJ_adaptive_mean, preco, runf, param, param_sub);

		for (i = 0; i<param->row; i++)
		for (j = 0; j<param->col; j++)
		{
			i1 = (int)(i + rp1[s]);
			j1 = (int)(j + rp1[s + f_dim]);
            
			if (i1 > 0 && i1 < param->row && j1 > 0 && j1 < param->col)
            {
				ind = i1 + param->row*j1;
				num_corrSurf = IJ_adaptive_mean[ind] - I_adaptive_mean[ind] * J_adaptive_mean[ind];
				dem_left = II_adaptive_mean[ind] - I_adaptive_mean[ind] * I_adaptive_mean[ind];
				dem_right = JJ_adaptive_mean[ind] - J_adaptive_mean[ind] * J_adaptive_mean[ind];
                dem_corrSurf = sqrt(dem_left*dem_right);
                
                if( dem_corrSurf > 1e-10)
                    fout[i + param->row*(j + param->col*s)] = exp(-(1-num_corrSurf/dem_corrSurf)/0.5);
                else
                    fout[i + param->row*(j + param->col*s)] = 1;
                
			}
		}
	}

	for (i = 0; i<param->row; i++)
	for (j = 0; j<param->col; j++)
	{
		fnorm = 0;
		for (s = 0; s<f_dim; s++)
			fnorm += fout[i + param->row*(j + param->col*s)] * fout[i + param->row*(j + param->col*s)];
        
        if( sqrt(fnorm) > 1e-10 )
            for (s = 0; s<f_dim; s++)
                fVol[i + row_o*(j + col_o*s)] = fout[i + param->row*(j + param->col*s)] / sqrt(fnorm);
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
    
	free(preco->integral_img);
	free(preco->glob_I_sub);
	free(preco->glob_I_mean);
	free(preco->glob_II_mean);
	free(preco->glob_I_var);
	free(preco);

	free(runf->a_mean);
	free(runf->b_mean);
	free(runf->p_mean);
	free(runf->Ip_mean);
	free(runf->p_sub);
	free(runf->a_mean_large);
	free(runf->b_mean_large);
	free(runf);

	free(diff_rp);
}

static inline void imresize_nearest(double *img, double *img_out, PARAM *param1, PARAM *param2)
{
	int i, j, x, y;
	double tx, ty;
	int row1 = param1->row, col1 = param1->col;
	int row2 = param2->row, col2 = param2->col;

	tx = (double)col1 / (double)col2;
	ty = (double)row1 / (double)row2;

	double *roop_img_out = img_out;
	for (j = 0; j<col2; j++)
	for (i = 0; i<row2; i++)
	{
		x = (int)(tx*j);
		y = (int)(ty*i);

		double *roop_img = img + y + row1*x;
		*roop_img_out = *roop_img;
		roop_img_out++;
	}
}

static inline void imresize_bilinear(double *img, double *img_out, PARAM *param1, PARAM *param2)
{
	int i, j, x, y;
	int ind;
	double tx, ty;
	double x_diff, y_diff, inv_x_diff, inv_y_diff;
	int row1 = param1->row, col1 = param1->col;
	int row2 = param2->row, col2 = param2->col;

	tx = (double)col1 / (double)col2;
	ty = (double)row1 / (double)row2;

	double *roop_img_out = img_out;
	for (j = 0; j<col2; j++)
	for (i = 0; i<row2; i++)
	{
		x = (int)(tx*j);
		y = (int)(ty*i);

		x_diff = (double)(tx*j) - x;
		y_diff = (double)(tx*i) - y;
		inv_x_diff = 1 - x_diff;
		inv_y_diff = 1 - y_diff;

		int col_bound = col1 - 2;
		int row_bound = row1 - 2;
		if (x>col_bound) x = col_bound;
		if (y>row_bound) y = row_bound;

		ind = y + row1*x;
		double *roop_img1 = img + ind;
		double *roop_img2 = img + ind + row1;
		double *roop_img3 = img + ind + 1;
		double *roop_img4 = img + ind + 1 + row1;

		*roop_img_out = *roop_img1 * inv_x_diff*inv_y_diff + *roop_img2 * x_diff*inv_y_diff
			+ *roop_img3 * inv_x_diff*y_diff +*roop_img4 * y_diff*x_diff;
		roop_img_out++;
	}
}

static inline void boxfilter(double *img, double *img_out, double *integral_img, PARAM *param_sub)
{
	int i,j;
	int init;
    int ind, ind1, ind2, ind3, ind4;
	int row = param_sub->row, col = param_sub->col, r = param_sub->r, len = param_sub->len;
	double win_size = (2 * r + 1)*(2 * r + 1);

    memset(integral_img, 0, len*sizeof(double));

	for( i=0;i<row;i++ )
		for( j=0;j<col;j++ )
        {
            ind = i + row*j;
            if( i==0 && j==0 )
                integral_img[ind] = img[ind];
            else if( i==0 && j!=0 )
                integral_img[ind] = img[ind] + integral_img[ind-row];
            else if( i!=0 && j==0 )
                integral_img[ind] = img[ind] + integral_img[ind-1];
            else
                integral_img[ind] = img[ind] + integral_img[ind-1] + integral_img[ind-row] - integral_img[ind-1-row];
        }
    
	for( i=r+1; i<row-r; i++ )
        for( j=r+1; j<col-r; j++ )
        {
            ind = i + row*j;
            ind1 = i+r + row*(j+r);
            ind2 = i-r-1 + row*(j+r);
            ind3 = i+r + row*(j-r-1);
            ind4 = i-r-1 + row*(j-r-1);
            
            img_out[ind] = (integral_img[ind1] - integral_img[ind2] - integral_img[ind3]
                    + integral_img[ind4]) / win_size;
        }
}

static inline void guidedfilter_precompute(double *I, PRECOMPUTE *preco, PARAM *param, PARAM *param_sub)
{
	int row = param_sub->row, col = param_sub->col, len = param_sub->len;

	memset(preco->glob_I_mean, 0, param_sub->len*sizeof(double));

	imresize_nearest(I, preco->glob_I_sub, param, param_sub);

	double *roop_glob_II_mean = preco->glob_II_mean;
	double *roop_glob_I_sub = preco->glob_I_sub;
	for (int i = 0; i < len; i++)
	{
		*roop_glob_II_mean = (*roop_glob_I_sub)*(*roop_glob_I_sub);
		roop_glob_II_mean++; 
		roop_glob_I_sub++;
	}

	boxfilter(preco->glob_I_sub, preco->glob_I_mean, preco->integral_img, param_sub);
	boxfilter(preco->glob_II_mean, preco->glob_II_mean, preco->integral_img, param_sub);

	double *roop_glob_I_var = preco->glob_I_var;
	double *glob_II_mean = preco->glob_II_mean;
	double *glob_I_mean = preco->glob_I_mean;
	for (int i = 0; i < len; i++)
	{
		*roop_glob_I_var = *glob_II_mean - (*glob_I_mean)*(*glob_I_mean) + param->epsil;
		roop_glob_I_var++; 
		glob_II_mean++; 
		glob_I_mean++;
	}
}

static inline void guidedfilter_runfilter(double *I, double *p, double *I_out, PRECOMPUTE *preco, RUNFILTER *runf, PARAM *param, PARAM *param_sub)
{
	memset(runf->a_mean, 0, param_sub->len*sizeof(double));
	memset(runf->b_mean, 0, param_sub->len*sizeof(double));
	memset(runf->p_mean, 0, param_sub->len*sizeof(double));
	memset(runf->Ip_mean, 0, param_sub->len*sizeof(double));

	imresize_nearest(p, runf->p_sub, param, param_sub);

	double *roop_Ip_mean = runf->Ip_mean;
	double *roop_glob_I_sub = preco->glob_I_sub;
	double *roop_p_sub = runf->p_sub;

	for (int i = 0; i < param_sub->len; i++)
	{
		*roop_Ip_mean = (*roop_glob_I_sub) * (*roop_p_sub);
		roop_Ip_mean++; 
		roop_glob_I_sub++; 
		roop_p_sub++;
	}

	boxfilter(runf->p_sub, runf->p_mean, preco->integral_img, param_sub);
	boxfilter(runf->Ip_mean, runf->Ip_mean, preco->integral_img, param_sub);

	double *roop_a_mean = runf->a_mean;
	double *roop_b_mean = runf->b_mean;
	roop_Ip_mean = runf->Ip_mean;
	double *roop_glob_I_mean = preco->glob_I_mean;
	double *p_mean = runf->p_mean;
	double *roop_glob_I_var = preco->glob_I_var;
	double *roop_p_mean = runf->p_mean;

	for (int i = 0; i < param_sub->len; i++)
	{
		*roop_a_mean = (*roop_Ip_mean - (*roop_glob_I_mean)*(*p_mean)) / *roop_glob_I_var;
		*roop_b_mean = *roop_p_mean - (*roop_a_mean)*(*roop_glob_I_mean);
		roop_a_mean++; 
		roop_b_mean++;
		roop_Ip_mean++; 
		roop_glob_I_mean++; 
		p_mean++;
		roop_glob_I_var++; 
		roop_p_mean++;
	}

	boxfilter(runf->a_mean, runf->a_mean, preco->integral_img, param_sub);
	boxfilter(runf->b_mean, runf->b_mean, preco->integral_img, param_sub);

	memset(runf->a_mean_large, 0, param->len*sizeof(double));
	memset(runf->b_mean_large, 0, param->len*sizeof(double));

	imresize_bilinear(runf->a_mean, runf->a_mean_large, param_sub, param);
	imresize_bilinear(runf->b_mean, runf->b_mean_large, param_sub, param);

	double *roop_I_out = I_out;
	double *roop_I = I;
	double *roop_a_mean_large = runf->a_mean_large;
	double *roop_b_mean_large = runf->b_mean_large;

	for (int i = 0; i < param->len; i++)
	{
		*roop_I_out = (*roop_a_mean_large) * (*roop_I) + (*roop_b_mean_large);
		roop_I_out++;
		roop_a_mean_large++;
		roop_I++;
		roop_b_mean_large++;
	}
}