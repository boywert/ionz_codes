/**
 * @file   ionz_main.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 21:01:59 2014
 * 
 * @brief  Main program
 * 
 */

#include "ion.h"


/** 
 * Main program
 * 
 * @param argc 
 * @param argv 
 * 
 * @return 
 */
main(int argc, char **argv) {
  FILE  *inp,*outpp;
  int ii, jj, kk,ll,jk,mm,sfac;
  float  Nnion,vaa,epsilon,box;
  int  temp,flag_sup;
  float dr,r_min,r_max;
  char file1[300],file2[300],num[50];
  int Nnion,N1,N2,N3;
  float *nion,xh1,dump,mass1,mass2;  
  int nhalo,output_flag,in_flag;
  float Radii,his_z,zval;
  double robar,robarhalo,vfac,*vion,*roion;
  float **rra,**vva,**halo,**data,*dummy,junk1=1.0,junk2=0.0;
  float *Radii_list;
  int n_radii;
  int *NjobsperTask;
  int *JobsTask;
  double t_start, t_stop;
  float *buffer, *buffer_final;
#ifdef CHUNKTRANSFER
  int mpi_buffer=1000000;
  int cur_len;
#endif
  char densfilename[2000], sourcefilename[2000],z_out[1000];
  fftw_real ***nh, ***ngamma, ****nxion;
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
  NTask = 1;
  ThisTask = 0;
#endif //PARALLEL

  if(argc != 4) {
    printf("Need 3 inputs\n");
    exit(0);
  }
  sprintf(densfilename,"%s",argv[1]);
  sprintf(sourcefilename,"%s",argv[2]);
  sprintf(z_out,"%s",argv[3]);
  if(ThisTask == 0) {
    system("date");
    printf("Start semi-numerical reionization process\n");
  }
  pi=4.0*atan(1.0);

  read_params("input.ionz");

  Nnion = input_param.Nnion;
  nion=(float*)calloc(Nnion,sizeof(float));
  for(ii=0;ii<Nnion;ii++) {
    nion[ii] = input_param.nion[ii];
  }
  vion=(double*)calloc(Nnion,sizeof(double));
  roion=(double*)calloc(Nnion,sizeof(double));
  N1 = input_param.N1;
  N2 = input_param.N2;
  N3 = input_param.N3;
  vomegam = input_param.omegam;
  vomegalam = input_param.omegalam;
  vomegab = input_param.omegab;
  
  /// Allocating memory to different arrays
  Setting_Up_Memory_For_ionz(Nnion, N1,N2,N3,nh, ngamma, nxion);

  t_start =Get_Current_time();

  /// Allocate buffer to store 3D array
  buffer = malloc(sizeof(float)*N1*N2*N3);

  /// Use Task:0 to read density
  if(ThisTask == 0) {
    read_density(densfilename,buffer,&robar,N1,N2,N3,vomegam,vomegab);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&robar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  unpack_3d_array_mpi_transfer(buffer,nh,N1,N2,N3);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(ThisTask == 0)  {
  read_sources(sourcefilename,buffer,&robarhalo,N1,N2,N3);  
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robarhalo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  unpack_3d_array_mpi_transfer(buffer,ngamma,N1,N2,N3);
  free(buffer);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  t_stop = Get_Current_time();

  /* Sanity MPI check */
  if(ThisTask == 1)
    printf("N1=%d N2=%d N3=%d\n",N1,N2,N3);
  if(ThisTask == 0)
    printf("reading in data %lf s\n",t_stop-t_start);

  //calculating max and min radius for smoothing in units of grid size
  r_min=1.;
  r_max=pow((1.*N1*N2*N3),(1./3.))/2.;

  Radii_list = malloc(sizeof(float)*constants.max_Nradii); 
  NjobsperTask = malloc(sizeof(float)*NTask);
  n_radii = make_radii_list(Radii_list,r_min,r_max);
  for(jj=0;jj<NTask;jj++) {
    NjobsperTask[jj] = n_radii/NTask;
    if(jj < n_radii%NTask)
      NjobsperTask[jj]++;
    if(jj == ThisTask) {
      JobsTask = malloc(sizeof(int)*NjobsperTask[jj]);
      for(ii=0;ii<NjobsperTask[jj];ii++) {
	JobsTask[ii] = ii*NTask+ThisTask;
      }
    }
  }
  // The max smoothing radius here is set as half of the diagonal of the box
  // This can be changed, one can choose a redshift dependent function instead
  // Or one can choose a model for redshift evolution of the mean free path of the UV photons
  // We are showing the most simple case here
  
  // printf("robar=%e  robarhalo=%e ratio= %e\n",robar,robarhalo,robar/robarhalo);

  //subgrid re-ionization

  for(jk=0;jk<Nnion;jk++) {
    //calculating avg. ionization frction
    vion[jk]=0.0;
    roion[jk]=0.0;

    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++) {
	  if(nh[ii][jj][kk]>nion[jk]*ngamma[ii][jj][kk]) {
	    nxion[jk][ii][jj][kk]=nion[jk]*ngamma[ii][jj][kk]/nh[ii][jj][kk];
	  }    
	  else {
	    nxion[jk][ii][jj][kk]=1.;
	  }
	  vion[jk]+=nxion[jk][ii][jj][kk];
	  roion[jk]+=nxion[jk][ii][jj][kk]*nh[ii][jj][kk];	  
	}
    vion[jk]/=(1.*N1*N2*N3);
    roion[jk]/=(float)(robar*N1*N2*N3);
    if(ThisTask == 0)
      printf("Subgrid: obtained vol. avg. x_ion=%e mass avg. x_ion=%e\n",vion[jk],roion[jk]);
 }

  if(ThisTask == 0)
    printf("Start semi-numerical reionization process\n");



#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_start = Get_Current_time();
  for(ii=0;ii<NjobsperTask[ThisTask];ii++) {
    reionization(Radii_list[JobsTask[ii]], nh, ngamma, nxion, nion, Nnion, N1, N2, N3 );    
  }  
  fftw_free(ngamma);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish reionizing process %lf s\n",t_stop-t_start);
  t_start = Get_Current_time();

#ifdef PARALLEL
  pack_4d_array_mpi_transfer(nxion,buffer,Nnion, N1, N2, N3);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_stop = Get_Current_time();
  if(ThisTask == 0) {
    printf("Finish packing data %lf s\n",t_stop-t_start);
#ifdef PARALLEL
    buffer_final = malloc(sizeof(float)*Nnion*N1*N2*N3);
#endif
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  t_start = Get_Current_time();
  ii = 0;
  while (ii*mpi_buffer < Nnion*N1*N2*N3) {
    cur_len = min(mpi_buffer, Nnion*N1*N2*N3-ii*mpi_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&buffer[ii*mpi_buffer],&buffer_final[ii*mpi_buffer],cur_len,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
    ii++;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish finding max:split %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#else
  t_start = Get_Current_time();
  MPI_Reduce(buffer, buffer_final, Nnion*N1*N2*N3, MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish finding max:whole %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif // PARALLEL

  if(ThisTask == 0) {
    for(jk=0;jk<Nnion;jk++) {
      //calculating avg. ionization frction
      vion[jk]=0.0;
      roion[jk]=0.0;
      // Defining the ionization map output file name
      // This is based on the value of nion assigned to it
      sprintf(file2,"xHI_map_%s_%4.2f",z_out,nion[jk]);
      printf("Saving %s\n",file2);
      // Writing the x_HI map in binary
      // In the begining 3 integers are written which defines the size
      // of the x_HI array
#ifdef PARALLEL
      
#endif

      inp=fopen(file2,"w");
      fwrite(&N1,sizeof(int),1,inp);
      fwrite(&N2,sizeof(int),1,inp);
      fwrite(&N3,sizeof(int),1,inp);
      
      
#ifdef PARALLEL
      fwrite(&buffer_final[jk*N1*N2*N3],sizeof(float),N1*N2*N3,inp);
#else
      fwrite(&buffer[jk*N1*N2*N3],sizeof(float),N1*N2*N3,inp);	
#endif
      for(ii=0;ii<N1*N2*N3;ii++) {
	    xh1=(1.-nxion[jk][ii][jj][kk]);
	    xh1=(xh1 >0.0)? xh1: 0.0;
	    vion[jk]+=xh1;
	    nxion[jk][ii][jj][kk]=xh1; // store x_HI instead of x_ion
	    nhs[ii][jj][kk]=xh1*nh[ii][jj][kk]; // ro_HI on grid
	    roion[jk]+=nhs[ii][jj][kk];	
      }
	    
      fclose(inp);
      roion[jk]/=(N1*N2*N3); // mean HI density in grid units
      
      vion[jk]/=(1.*N1*N2*N3); // volume avg xHI
      roion[jk]/=(float)robar; // divide by H density to get mass avg. xHI
      printf("nion = %f obtained vol. avg. x_HI=%e mass avg. x_HI=%e\n",nion[jk],vion[jk],roion[jk]);
      system("date");
    }
  }
  MPI_Finalize();
  return 0;
} /* END OF MAIN */

