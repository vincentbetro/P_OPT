#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "geometry.h"
#include "Point.h"
#include "Vector.h"
#include "P_OPT.h"
#include "POLYMESH.h"
#include "Util.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

void smooth(geometry *geom, char sname[], int nsmoo, int mvbnd, int mvint, int sflag, int eflag, double threshold, int blend_exponent)
{
  const int bdim = 80;
  char buff[bdim];
  char *ptr;

  POLYMESH *mesh = new POLYMESH();

  if (num_procs > 1)
  {
    sprintf(buff,"%s",sname);
    ptr = strstr(buff,".sg");
    if (ptr == NULL)
    {
      fprintf(out_f,"\nSG suffix <.sg> not found in file name!");
      fflush(out_f);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
    } else
    {
      //reset cursor to overwrite extension with new extension and proc num
      *ptr = '\0';
    }
    sprintf(sname,"%s_%d.sg",buff,my_rank);
  }

  // read mesh data
  mesh->read_mesh(sname);
  
  mesh->create_node_lists(0);

  // examine current mesh quality
  mesh->mesh_stats();

  // call optimization-based smoother
  mesh->OPT_move(geom,nsmoo,mvbnd,mvint,sflag,eflag,threshold,blend_exponent);

  // overwrite mesh data
  mesh->write_mesh(sname);

  // examine current mesh quality
  mesh->mesh_stats();

  //double **func;
  //char (*func_name)[80];
  //int num_func = 0;
  //mesh->Fieldview(1,sname,num_func,func_name,func);
  //mesh->Ensight(1,sname,num_func,func_name,func);

  delete mesh;

}
