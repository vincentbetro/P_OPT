#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "geometry.h"
#include "P_OPT.h"
#include "Point.h"
#include "Vector.h"
#include "journal.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

//global variables
int my_rank; /*rank of process*/
int num_procs; /*number of processes*/
FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

int main(int argcs, char* pArgs[])
{
  int source; //rank of sender
  int tag = 0; //tag for messages
  int i, nf, nsmoo;
  int mvbnd, mvint, eflag, sflag, blend_exponent;
  double threshold;
  const int bdim = 132; //buffer dim
  char** fnames; //file name storage
  char buff[bdim]; //buffer
  char sname[bdim]; //mesh name storage
  char param_file[bdim]; // parameter file name
  time_t tm;
  clock_t time_0, time_1;
  char *t_char;
 
#ifdef PARALLEL
  int dest; //rank of receiver
  MPI_Status status; //return status for receive
  //Start up MPI
  MPI_Init(&argcs, &pArgs);
   
  //Find out process rank of current instance
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
  //Find out number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#else
  my_rank = 0;
  num_procs = 1;
#endif
  
  in_f = stdin;
  out_f = stdout;

  if (my_rank == 0)
  {
    //create journal file
    if ((jou_f=fopen("P_OPT.jou","w")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file journal file");
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    }
  }
  
  //check for standard input
  if (--argcs < 1)
  {
    if (my_rank == 0)
    {
      fprintf(stdout,"\nNo input file specified!");
      fprintf(stdout,"\nUsing standard input!");
      fflush(stdout);
    }
  } else
  {
    if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file <%s>\n",pArgs[argcs]);
      if (my_rank == 0) fprintf(jou_f,"\nCouldn't open file <%s>\n",pArgs[argcs]);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    }
  }

  if (in_f != stdin)
  {
    sprintf(buff,"P_OPT.out");

    if (my_rank == 0)
    {
      if ((out_f=fopen(buff,"w")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file output file %s",buff);
        fclose(jou_f);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
    } else
    {
      if ((out_f=fopen(buff,"wa")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file output file %s",buff);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
    }
  }
  if (my_rank != 0)
    fclose(in_f);

  // ALL PROCESSES NOW HAVE THE OUTPUT FILE AVAILABLE FOR PRINTING
  // ALL NORMAL PRINTING IS DONE BY PROCESS 0, OTHERS PRINT DEBUG INFORMATION
  // ALL INPUT IS THROUGH PROCESS 0

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // output the beginning wall time
  time(&tm);
  t_char = ctime(&tm);
  time_0 = clock();
  if (my_rank == 0)
  {
    fprintf(out_f,"\nP_OPT run started at %s",t_char);

    //print program info 
    fprintf(out_f,"\n====================================================================");
    fprintf(out_f,"\n COPYRIGHT 2003-2012 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA     ");
    fprintf(out_f,"\n                                                                    ");
    fprintf(out_f,"\n                    RIGHTS IN DATA                                  ");
    fprintf(out_f,"\n                                                                    ");
    fprintf(out_f,"\n THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT ");
    fprintf(out_f,"\n    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO       ");
    fprintf(out_f,"\n         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL      ");
    fprintf(out_f,"\n              RIGHTS IN DATA CONTRACT CLAUSES.                      ");
    fprintf(out_f,"\n       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY   ");
    fprintf(out_f,"\n              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA            ");
    fprintf(out_f,"\n                                                                    ");
    fprintf(out_f,"\n Parallel Optimization-based smoothing (P_OPT)                      ");
    fprintf(out_f,"\n NOTE: This data includes the UT SimCenter at Chattanooga P_OPT     ");
    fprintf(out_f,"\n code, which was developed under private non-government funding.    ");
    fprintf(out_f,"\n This software is submitted with limited rights to use, reproduce,  ");
    fprintf(out_f,"\n and disclose this data for Government Purposes only.               ");
    fprintf(out_f,"\n Requests for access to the software for non-governmental purposes  ");
    fprintf(out_f,"\n should be referrred to                                             "); 
    fprintf(out_f,"\n                                                                    ");
    fprintf(out_f,"\n    Dr. Steve Karman                                                "); 
    fprintf(out_f,"\n    Steve-Karman@utc.edu                                            "); 
    fprintf(out_f,"\n    423-425-5492  or  423-425-5470                                  "); 
    fprintf(out_f,"\n                                                                    ");
    fprintf(out_f,"\n    University of Tennessee SimCenter at Chattanooga                "); 
    fprintf(out_f,"\n    701 East M. L. King Boulevard                                   "); 
    fprintf(out_f,"\n    Chattanooga, TN 37403                                           "); 
    fprintf(out_f,"\n====================================================================\n");
    fflush(out_f);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (my_rank == 0)
  {
    // obtain number of geometry files
    journal(in_f, out_f, jou_f, "#Number of geometry files ->",nf);

    if (nf <= 0)
    {
      fprintf(out_f,"\nNumber of geometry files must be > 0\n");
      fflush(out_f);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    }
	  
#ifdef PARALLEL
    // inform slave processes
    tag = 0;
    for (dest = 1; dest < num_procs; dest++)
      MPI_Send(&nf, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
#endif

    // obtain geometry file names
    fnames = (char**)malloc(nf*sizeof(char*));
    
    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));
      sprintf(buff,"#File %d >",i+1);
      journal(in_f, out_f, jou_f, buff, fnames[i]);

#ifdef PARALLEL
      // inform slave processes
      tag = i+1;
      for (dest = 1; dest < num_procs; dest++)
        MPI_Send(fnames[i], strlen(fnames[i])+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
#endif

    }
  } else
  {
    // Obtain number of geometry files from master
    source = 0;
    tag = 0;
#ifdef PARALLEL
    MPI_Recv(&nf, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    // Obtain geometry file names from master
    fnames = (char**)malloc(nf*sizeof(char*));

    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));

      tag = i+1;
      MPI_Recv(fnames[i], bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
    }
#endif
  }

  // read geometry from file list
  int *facet = 0;
  int *mat = 0;
  double *gvert = 0;
  int ngn = 0;
  int nfacets = 0;
  read_geometry(nf, fnames, ngn, nfacets, gvert, facet, mat, out_f);
  
  sprintf(param_file,"%s","geometry.params");
  geometry* geom = new geometry();
  geom->process_geom(ngn, nfacets, gvert, facet, mat, 0, 0, 0, param_file, out_f);

  // output the geometry prep wall time for master
  time(&tm);
  t_char = ctime(&tm);
  time_1 = clock();
  if (my_rank == 0)
  {
    fprintf(out_f,"\nTotal geometry prep time in seconds = %14.7e\n",(float)(time_1-time_0)/CLOCKS_PER_SEC);
    fflush(out_f);
  }
  
  //clean up fnames
  for (i=0; i < nf; i++)
    free(fnames[i]);
  free(fnames);
  
  if (my_rank == 0)
  {
    //read in grid file name
    journal(in_f, out_f, jou_f, "#Enter physical grid file name (serial version) ->",sname);
  
#ifdef PARALLEL
    // inform slave processes of file name, using strlen for count
    tag = my_rank;
    for (dest = 1; dest < num_procs; dest++)
      MPI_Send(sname, strlen(sname)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);	
#endif
  } else
  {
    source = 0; //reset to master
    tag = 0;
#ifdef PARALLEL
    MPI_Recv(sname, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
#endif
  }
    
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  double pi=4*atan(1.0);
  if (my_rank == 0)
  {
    journal(in_f, out_f, jou_f, "#Enter # of smoothing sweeps > ",nsmoo);
    journal(in_f, out_f, jou_f, "#Enter boundary movement flag [0 1] > ",mvbnd);
    journal(in_f, out_f, jou_f, "#Enter interior movement flag [<0, 0, >0] > ",mvint);
    journal(in_f, out_f, jou_f, "#Enter threshold parameter [0.0-1.0] > ", threshold);
    journal(in_f, out_f, jou_f, "#Enter surface cost flag [0-no, 1-yes, 2-mixed] > ", sflag);
    journal(in_f, out_f, jou_f, "#Enter expanded stencil flag [0-no, 1-yes, 2-mixed] > ", eflag);
    journal(in_f, out_f, jou_f, "#Cost function blend exponent [ >= 0] > ", blend_exponent);

    if (in_f != stdin) fclose(in_f);

    fclose(jou_f);
  }

#ifdef PARALLEL
  MPI_Bcast(&nsmoo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mvbnd, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mvint, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&threshold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&eflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&blend_exponent, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  if (my_rank == 0)
  {
    fprintf(out_f,"\n\nInput parameters:");
    fprintf(out_f,"\nNumber of smoothing sweeps = %d",nsmoo);
    fprintf(out_f,"\nBoundary movement flag = %d",mvbnd);
    fprintf(out_f,"\nInterior movement flag = %d",mvint);
    fprintf(out_f,"\nThreshold parameter = %g",threshold);
    fprintf(out_f,"\nSurface cost flag = %d",sflag);
    fprintf(out_f,"\nExpanded stencil flag = %d",eflag);
    fprintf(out_f,"\nCost function blend exponent = %d",blend_exponent);
    fprintf(out_f,"\n");
    fflush(out_f);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  smooth(geom,sname,nsmoo,mvbnd,mvint,sflag,eflag,threshold,blend_exponent);
 
  //perform end time check
  time(&tm);
  t_char = ctime(&tm);
  time_1 = clock();
  if (my_rank == 0)
    fprintf(out_f,"\nTotal smoothing time in seconds = %14.7e\n",(float)(time_1-time_0)/CLOCKS_PER_SEC);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  delete geom;

  if (my_rank == 0)
    fclose(out_f);

#ifdef PARALLEL
  MPI_Finalize(); 
#endif
  
  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0)
    fprintf(out_f,"\nP_OPT run completed at %s",t_char);

  if (out_f != stdout) fclose(out_f);

  return(0);
}
