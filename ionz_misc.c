/**
 * @file   ionz_misc.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 23:55:03 2014
 * 
 * @brief  Misc functions
 * 
 * 
 */

#include "ion.h"

float cal_H(float z, float H0, float Om) {
  return H0*sqrt(Om*(1.+z)*(1.+z)*(1.+z) + (1.-Om));
}
double delta_t(float z_max, float z_min,float Om, float h) {
  int i,n = 100000;
  const float H0 = 100.0;
  const float Mpc2m = 3.08567758e22;
  double dt,dz,F1,F2,z1,z2;
  dt = 0.;
  dz = (z_max - z_min)/n;
  for(i=0;i<n;i++) {
    z1 = z_min + dz*i;
    z2 = z1+dz;
    F1 = 1./(1.+z1)/cal_H(z1,H0,Om);
    F2 = 1./(1.+z2)/cal_H(z2,H0,Om);
    dt += 0.5*(F1+F2)*dz;
  }
  return dt*(Mpc2m/h/1000.);
}


void chunk_float_mpi_bcast(float *buffer, int len, int mpi_buffer,int root) {
  int ii, cur_len;
  ii = 0;
  while(ii*mpi_buffer < len) {
    cur_len = min(mpi_buffer, len-ii*mpi_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&buffer[ii*mpi_buffer], cur_len, MPI_FLOAT, 0, MPI_COMM_WORLD);
    ii++;
  }   
}

/** 
 * Free 3D fftw_real arrays
 * 
 * @param f fftw_real*** to free 
 * @param N1 1st dimension len 
 * @param N2 2nd dimension len
 * @param N3 3rd dimension len
 */
void free_fftw_real_3d(fftw_real ***f, int N1, int N2, int N3) {
  int ii;
  free(&f[0][0][0]);
  for(ii=0;ii<N1;ii++) 
    free(f[ii]);
  free(f);
}

/** 
 * Allocate 3D fftw_real arrays
 * 
 * @param N1 1st dimension len 
 * @param N2 2nd dimension len
 * @param N3 3rd dimension len
 * 
 * @return 3D fftw_real array pointer
 */
fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3) {
  int ii,jj;
  fftw_real ***phia,*phi;
 
  phia=(fftw_real ***)malloc(N1 *  sizeof(fftw_real **));
  phi = (fftw_real *) calloc(N1 * N2 * N3,sizeof(fftw_real));

  for(ii=0;ii<N1;ii++) {
      phia[ii]=(fftw_real **)malloc(N2 *  sizeof(fftw_real*));
      for(jj=0;jj<N2;jj++)
	phia[ii][jj]= phi+ (N2*N3)*ii + N3*jj;
  }
  return(phia);
}

/** 
 * Get current time
 * 
 * 
 * @return  current time
 */
double Get_Current_time() {
#ifdef PARALLEL
  double time;
  time = MPI_Wtime();
  return time;
#else
  time_t timer;
  timer = time(NULL);
  return (double) timer;
#endif
}

/** 
 * Generate a list of radii to compute semi-numerical reionization 
 * 
 * @param radii_p Pointer to radii list (pre-allocated)
 * @param r_min Minimum radius
 * @param r_max Maximum radius
 * @param dr_inc Increase in radius (0.0 < dr_inc <= 1.0; 0.1 by default)
 * @param max_dr Maximum increase in radius (2.0 by default)
 * @return total number of radii 
 */
int make_radii_list(float *radii_p, float r_min, float r_max, float dr_inc, float max_dr) {
  float r,dr;
  int i;

  /// Test if nRadii <= max_Nradii
  i = 0;
  r = r_min;
  while (r < r_max) {
    dr = min(r*dr_inc,max_dr);
    r += dr;
    i++;
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(mympi.ThisTask == 0) 
    if(i > constvars.max_Nradii) {
      printf("[%s][%s] ***Error: File:%s Line:%d\n",__DATE__,__TIME__,__FILE__,__LINE__);
      printf("The number of smoothing radii is greater than %d\nTerminate\n",constvars.max_Nradii);
      exit(1);
    }
  r = r_min;
  i = 0;
  /// Store radii list
  while (r < r_max) {
    radii_p[i] = r;
    dr = min(r*dr_inc,max_dr);
    r += dr;
    i++;
  }  
  radii_p = realloc(radii_p,sizeof(float)*i);
  return i;
}
