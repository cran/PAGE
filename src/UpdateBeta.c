
# include <R.h>
# include <Rmath.h>
void UpdateOnebeta(double *Y, double *X, double *B, int *dims, double *C,
double *W, double *lambda, int *idx);

void UpdateBeta(double *Y, double *X, double *B, int *dims, double *C,
double *W, double *lambda, double *updateB, int *iteration)
{
    int i=1, j;
    int M = dims[1];
    int P = dims[2];
    int idx[2];
    int iter_n = 0;
    double maxDiff = 1;
    int iter_count = 0;
    
    while(maxDiff>1e-4 && iter_n<1000)
    {
        maxDiff = 0;
        
        for(i=0;i<M;i++)
        {
          for(j=0;j<P;j++)
          {
            idx[0] = i;
            idx[1] = j;

            double old_B = B[j*M+i];
            double tempDiff;

            UpdateOnebeta(Y, X, B,dims,C,W,lambda,idx);
	    iter_count++;
            tempDiff = sign(old_B-B[j*M+i])*(old_B-B[j*M+i]);
            maxDiff = (maxDiff>tempDiff)? maxDiff:tempDiff;
          }
        }
      iter_n++;
    }
    
    for(i=0;i<M;i++)
      for(j=0;j<P;j++)
        updateB[j*M+i] = B[j*M+i];
    *iteration = iter_count;
}



void UpdateOnebeta(double *Y, double *X, double *B, int *dims, double *C,
double *W, double *lambda, int *idx)
{
    int i, j, k;
    int N = dims[0];
    int M = dims[1];
    int P = dims[2];
    int m = idx[0];
    int p = idx[1];
    double *residual = (double *)malloc(N*M*sizeof(double));
    double *presidual = (double *)malloc(M*sizeof(double));
    double sign_factor1,sign_factor2,sign_factor, updated_beta;
    
    for(i=0;i<N;i++)
        for(j=0;j<M;j++)
        {
          residual[i*M+j]=0;
          for(k=0;k<P;k++)
              residual[i*M+j]+=X[i*P+k]*B[k*M+j];
	  residual[i*M+j]=Y[i*M+j]-residual[i*M+j];
        }
    for(i=0;i<M;i++)
    {
	     presidual[i]=0;
	     for(k=0;k<N;k++)
		      presidual[i]+=X[k*P+p]*residual[k*M+i];
    }

    for(i=0,sign_factor1=0;i<M;i++)
    	 sign_factor1+=C[i*M+m]*presidual[i];

    for(i=0,sign_factor2=0;i<N;i++)
	     sign_factor2+=X[i*P+p]*X[i*P+p];
    sign_factor2 = C[m*M+m]*sign_factor2;

    sign_factor = sign_factor1/sign_factor2+B[p*M+m];

    updated_beta = sign(sign_factor)*sign_factor-(*lambda)*W[p*M+m]/(2*sign_factor2);
    B[p*M+m] = sign(sign_factor)*updated_beta*(updated_beta>0);
    
    free(residual);
    free(presidual);
}
