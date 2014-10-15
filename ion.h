/**
 * @file   ion.h
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 20:49:08 2014
 * 
 * @brief  
 * Define all functions
 * 
 */
#ifndef ION_H_

#include <srfftw.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif


#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))

#define debug_checkpoint() (printf("CHECKPOINT: Fuction:%s Line:%d File:%s\n",__func__,__LINE__,__FILE__))

#define PREFIX "xHI_map_"

/// Global variables
typedef struct CONSTANTS {
  /// PI
  float pi;
  /// Maximum number of radii
  int max_Nradii;
  /// Increase in radius (factor of radius)
  float dr_inc;
  /// Maximum increase in radius
  float max_dr;
} struct_const;
extern struct_const constvars;


/// MPI variables
struct myMPI {
  /// Total MPI rank
  int NTask;
  /// This MPI node number
  int ThisTask;
} mympi;


/// Input Parameters
struct params 
{
  /// Number of Nions used
  int Nnion;
  /// Ngamma/Nh 
  float nion[100];
  /// expansion factor
  float a_expansion;
  /// redshift
  float z;
  /// Hubble_h constant 
  float Hubble_h;
  /// Omega_matter
  float omegam;
  /// Omega_Lambda
  float omegalam;
  /// Omega_baryon
  float omegab;
  /// 1st dimension grid
  int N1;
  /// 2nd dimension grid
  int N2;
  /// 3rd dimension grid
  int N3;
  /// Boxsize (Mpc/h)
  float boxsize;
  /// Gridsize (Mpc/h)
  float gridsize;
  /// Density file
  char densityfile[2000];
  /// sources file
  char sourcesfile[2000];
  /// current redshift
  char cur_z[100];
  /// previous redshift
  char prev_z[100];
  /// Output folder
  char outputdir[2000];
} input_param;

// in read_param.c
extern void read_nion(char *filename);
extern void read_params(char *filename);

/* in ionz_misc.c */
extern int make_radii_list(float *radii_p, float r_min, float r_max, float dr_inc, float max_dr);
extern double Get_Current_time();

/* in ionz_io.c */
extern void pack_3d_array_mpi_transfer(fftw_real ***input, float *output, int N1, int N2, int N3);
extern void unpack_3d_array_mpi_transfer(float *input, fftw_real ***output, int N1, int N2, int N3);
extern void pack_4d_array_mpi_transfer(fftw_real ****input, float *output, int Nnion, int N1, int N2, int N3);
extern void unpack_4d_array_mpi_transfer(float *input, fftw_real ****output, int Nnion,int N1, int N2, int N3);
extern void read_density(char *filename, float *buffer_3d, double *robar_p, int N1, int N2, int N3, float vomegam, float vomegab);
extern void read_xfrac(char *dirname, float *buffer_4d, float *nion_list, int Nnion, int N1, int N2, int N3);
extern void read_sources(char *filename, float *buffer_3d, double *robarhalo_p, int N1, int N2, int N3);

/* in allotarrays.c */
extern fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3);
extern void free_fftw_real_3d(fftw_real ***f, int N1, int N2, int N3);
extern void bcast_input_params();

/* in ionz_funcs.c */
extern void Setting_Up_Memory_For_ionz(int Nnion, int N1, int N2, int N3);
extern void smooth(fftw_real ***ro_dum,float Radii,int N1,int N2, int N3);
extern void subgrid_reionization(fftw_real ***nh, fftw_real ***ngamma, fftw_real ****nxion, double robar,float *nion, int Nnion, int N1, int N2, int N3);
extern void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3);
extern void subgrid_reionization_with_xfrac(fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****xfrac_p, fftw_real ****nxion_p, double robar, float *nion_p, int Nnion, int N1, int N2, int N3);
extern void reionization_with_xfrac(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****xfrac_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3);
/// Global arrays


#define ION_H_
#endif










