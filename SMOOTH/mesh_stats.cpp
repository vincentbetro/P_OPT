#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "Point.h"
#include "Vector.h"
#include "smooth.h"
#include "Util.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

double POLYMESH::element_volume(int e, Point &cg)
{
  Point p0, p1, p2;
  double vol;
  int n, n0, n1, n2, f, i, j;

  cg = Point(0.0,0.0,0.0);
  n=0;
  for (j=0; j < element[e].face_list->max; j++)
  {
    f = abs(element[e].face_list->list[j]);
    for (i=0; i < face[f].node_list->max; i++)
    {
      cg += node[face[f].node_list->list[i]].vert;
      n++;
    }
  }
  cg /= MAX(1,n);

  vol = 0.0;
  for (j=0; j < element[e].face_list->max; j++)
  {
    f = element[e].face_list->list[j];
    n0 = face[abs(f)].node_list->list[0];
    p0 = node[n0].vert;
    n1 = face[abs(f)].node_list->list[1];
    for (i=2; i < face[abs(f)].node_list->max; i++)
    {
      n2 = face[abs(f)].node_list->list[i];
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      if (f > 0)
        vol += tetrahedral_volume(cg, p0, p1, p2);
      else
        vol -= tetrahedral_volume(cg, p0, p1, p2);
      n1 = n2;
    }
  }

  return(vol);
}

void POLYMESH::mesh_stats()
{
  int b, count, e, f, i, j, k, local, n, n0, n1, n2, n3, n4, n5, n6, n7;
  Vector v1, v2, v3;
  Point cg, p0, p1, p2, p3, p4, p5, p6, p7, pmin, pmax, qmin, qmax;
  double ds, vol, fmin, favg, fmax, dlocal, mag, javg, jmin, jmax;
  int gnn, gntet, gnpyr, gnpri, gnhex, gnply, gnt, gnq, gngon;
  int neg, neg_skw, pos_skw, pos;
  double wgt[3][3];

  double pi = 4.0*atan(1.0);
#ifdef PARALLEL
    struct 
    { 
      double val; 
      int proc; 
    } g_in, g_out; 
    int globali;
    double globald;
    MPI_Status single_status;

    MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (my_rank == 0) fprintf(out_f,"\n\nMesh Statistics:");

  local=0;
  for (n=0; n < nn; n++)
    if (node[n].me.proc == my_rank) local++;
  gnn = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnn,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of nodes = %i",gnn);

  local=0;
  for (e=0; e < nelems; e++)
    if (element_type(e) == 4 && element[e].me.proc == my_rank) local++;
  gntet = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gntet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of tetrahedra = %i",gntet);

  local=0;
  for (e=0; e < nelems; e++)
    if (element_type(e) == 5 && element[e].me.proc == my_rank) local++;
  gnpyr = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of pyramids = %i",gnpyr);

  local=0;
  for (e=0; e < nelems; e++)
    if (element_type(e) == 6 && element[e].me.proc == my_rank) local++;
  gnpri = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of prisms = %i",gnpri);

  local=0;
  for (e=0; e < nelems; e++)
    if (element_type(e) == 8 && element[e].me.proc == my_rank) local++;
  gnhex = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnhex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of hexahedra = %i",gnhex);

  local=0;
  for (e=0; e < nelems; e++)
    if (element_type(e) == 0 && element[e].me.proc == my_rank) local++;
  gnply = local;
#ifdef PARALLEL
    MPI_Allreduce(&local,&gnply,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nTotal number of polyhedra = %i",gnply);

  if (my_rank == 0) fprintf(out_f,"\nTotal number of boundaries = %d",nb);

  // boundary statistics

  for (b=0; b < nb; b++)
  {
    if (my_rank == 0) fprintf(out_f,"\n\nStatistics for boundary %d : %s",b+1,boundary[b].name);

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count=0;
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = abs(boundary[b].polygon_list->list[i]);
      for (j=0; j < face[f].node_list->max; j++)
      {
        n1 = face[f].node_list->list[j];
        if (j < face[f].node_list->max-1)
          n2 = face[f].node_list->list[j+1];
        else
          n2 = face[f].node_list->list[0];
        ds = distance(node[n1].vert,node[n2].vert);
        fmin = MIN(fmin,ds);
        favg += ds;
        fmax = MAX(fmax,ds);
        count++;
      }
    }

#ifdef PARALLEL
    dlocal = fmin;
    MPI_Allreduce(&dlocal,&fmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMinimum tangential spacing = %g",fmin);

#ifdef PARALLEL
    local = count;
    MPI_Allreduce(&local,&count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dlocal = favg;
    MPI_Allreduce(&dlocal,&favg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    favg /= MAX(1,count);
    if (my_rank == 0) fprintf(out_f,"\nAverage tangential spacing = %g",favg);

#ifdef PARALLEL
    dlocal = fmax;
    MPI_Allreduce(&dlocal,&fmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMaximum tangential spacing = %g",fmax);

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count=0;
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = boundary[b].polygon_list->list[i];
      Vector norm = Vector(0.0,0.0,0.0);
      for (j=2; j < face[abs(f)].node_list->max; j++)
      {
        n1 = face[abs(f)].node_list->list[0];
        n2 = face[abs(f)].node_list->list[j-1];
        n3 = face[abs(f)].node_list->list[j];
        v1 = Vector(node[n1].vert,node[n2].vert);
        v2 = Vector(node[n1].vert,node[n3].vert);
        if (f > 0)
          norm += v1 % v2;
        else
          norm += v2 % v1;
      }
      norm.normalize();

      for (j=0; j < node[n1].element->max; j++)
      {
        e = node[n1].element->list[j];
        if (!element[e].face_list->Is_In_List(-f)) continue;
        for (k=0; k < element[e].node_list->max; k++)
        {
          n = element[e].node_list->list[k];
          if (face[abs(f)].node_list->Is_In_List(n)) continue;
          v1 = Vector(node[n1].vert,node[n].vert);
          ds = v1*norm;
          fmin = MIN(fmin,ds);
          favg += ds;
          fmax = MAX(fmax,ds);
          count++;
        }
      }
    }

#ifdef PARALLEL
    dlocal = fmin;
    MPI_Allreduce(&dlocal,&fmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMinimum normal spacing = %g",fmin);

#ifdef PARALLEL
    local = count;
    MPI_Allreduce(&local,&count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dlocal = favg;
    MPI_Allreduce(&dlocal,&favg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    favg /= MAX(1,count);
    if (my_rank == 0) fprintf(out_f,"\nAverage normal spacing = %g",favg);

#ifdef PARALLEL
    dlocal = fmax;
    MPI_Allreduce(&dlocal,&fmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMaximum normal spacing = %g",fmax);

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count=0;
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = boundary[b].polygon_list->list[i];
      Vector norm = Vector(0.0,0.0,0.0);
      for (j=2; j < face[abs(f)].node_list->max; j++)
      {
        n1 = face[abs(f)].node_list->list[0];
        n2 = face[abs(f)].node_list->list[j-1];
        n3 = face[abs(f)].node_list->list[j];
        v1 = Vector(node[n1].vert,node[n2].vert);
        v2 = Vector(node[n1].vert,node[n3].vert);
        if (f > 0)
          norm += v1 % v2;
        else
          norm += v2 % v1;
      }
      mag = norm.magnitude();

      fmin = MIN(fmin,mag);
      favg += mag;
      fmax = MAX(fmax,mag);
      count++;
    }

#ifdef PARALLEL
    dlocal = fmin;
    MPI_Allreduce(&dlocal,&fmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMinimum area = %g",fmin);

#ifdef PARALLEL
    local = count;
    MPI_Allreduce(&local,&count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dlocal = favg;
    MPI_Allreduce(&dlocal,&favg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    favg /= MAX(1,count);
    if (my_rank == 0) fprintf(out_f,"\nAverage area = %g",favg);

#ifdef PARALLEL
    dlocal = fmax;
    MPI_Allreduce(&dlocal,&fmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMaximum area = %g",fmax);

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count=0;
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = boundary[b].polygon_list->list[i];
      for (j=0; j < face[abs(f)].node_list->max; j++)
      {
        n1 = face[abs(f)].node_list->list[j];
        if (j < face[abs(f)].node_list->max-1)
          n2 = face[abs(f)].node_list->list[j+1];
        else
          n2 = face[abs(f)].node_list->list[0];
        if (j > 0)
          n3 = face[abs(f)].node_list->list[j-1];
        else
          n3 = face[abs(f)].node_list->list[face[abs(f)].node_list->max-1];
        v1 = Vector(node[n1].vert,node[n2].vert);
        v2 = Vector(node[n1].vert,node[n3].vert);
        v1.normalize();
        v2.normalize();
        mag = v1*v2;
        mag = MAX(-1.0,MIN(1.0,mag));
        mag = acos(mag)*180.0/pi;

        fmin = MIN(fmin,mag);
        favg += mag;
        fmax = MAX(fmax,mag);
        count++;
      }
    }

#ifdef PARALLEL
    dlocal = fmin;
    MPI_Allreduce(&dlocal,&fmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMinimum angle (degrees) = %g",fmin);

#ifdef PARALLEL
    local = count;
    MPI_Allreduce(&local,&count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    dlocal = favg;
    MPI_Allreduce(&dlocal,&favg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    favg /= MAX(1,count);
    if (my_rank == 0) fprintf(out_f,"\nAverage angle (degrees) = %g",favg);

#ifdef PARALLEL
    dlocal = fmax;
    MPI_Allreduce(&dlocal,&fmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (my_rank == 0) fprintf(out_f,"\nMaximum angle (degrees) = %g",fmax);
  }

  // element statistics

  fmin = 1.0e20;
  favg = 0.0;
  fmax = -1.0e20;
  count = 0;
  for (e=0; e < nelems; e++)
  {
    if (element[e].me.proc != my_rank) continue;
    count++;

    vol = element_volume(e,cg);

    if (vol < fmin)
    {
      fmin = vol;
      pmin = cg;
    }
    favg += vol;
    if (vol > fmax)
    {
      fmax = vol;
      pmax = cg;
    }
  }
  favg /= MAX(1,count);

  if (my_rank == 0) fprintf(out_f,"\n");

  if (num_procs > 1)
  {
#ifdef PARALLEL
    g_in.val = fmin;
    g_in.proc = my_rank;
    MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (my_rank == 0)
    {
      if (g_out.proc != my_rank)
      {
        MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
      }
    } else
    {
      if (g_out.proc == my_rank)
      {
        MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
      }
    }

    favg *= count;
    MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
    MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
    favg = globald/globali;

    g_in.val = fmax;
    g_in.proc = my_rank;
    MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (my_rank == 0)
    {
      if (g_out.proc != my_rank)
      {
        MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
      }
    } else
    {
      if (g_out.proc == my_rank)
      {
        MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
      }
    }
#endif
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nMinimum element volume = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
    fprintf(out_f,"\nAverage element volume = %lg",favg);
    fprintf(out_f,"\nMaximum element volume = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
  }

  if (gntet > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 4) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      cg = (p0+p1+p2+p3)/4.0;

      mag = tetrahedral_aspect_ratio(p0,p1,p2,p3);

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral aspect ratio = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage tetrahedral aspect ratio = %lg",favg);
      fprintf(out_f,"\nMaximum tetrahedral aspect ratio = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnpyr > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 5) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      cg = (p0+p1+p2+p3+p4)/5.0;

      mag = pyramid_aspect_ratio(p0,p1,p2,p3,p4);

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid aspect ratio = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage pyramid aspect ratio = %lg",favg);
      fprintf(out_f,"\nMaximum pyramid aspect ratio = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnpri > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 6) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;
      cg = (p0+p1+p2+p3+p4+p5)/6.0;

      mag = prism_aspect_ratio(p0,p1,p2,p3,p4,p5);

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism aspect ratio = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage prism aspect ratio = %lg",favg);
      fprintf(out_f,"\nMaximum prism aspect ratio = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnhex > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 8) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      n6 = element[e].node_list->list[6];
      n7 = element[e].node_list->list[7];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;
      p6 = node[n6].vert;
      p7 = node[n7].vert;
      cg = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;

      mag = hexahedral_aspect_ratio(p0,p1,p2,p3,p4,p5,p6,p7);

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral aspect ratio = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage hexahedral aspect ratio = %lg",favg);
      fprintf(out_f,"\nMaximum hexahedral aspect ratio = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gntet > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 4) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      cg = (p0+p1+p2+p3)/4.0;

      mag = tetrahedral_condition_number(p0,p1,p2,p3,1);

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral weighted condition number = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage tetrahedral weighted condition number = %lg",favg);
      fprintf(out_f,"\nMaximum tetrahedral weighted condition number = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnpyr > 0)
  {
    wgt[0][0] = 1.0;
    wgt[1][0] = 0.0;
    wgt[2][0] = 0.0;
    wgt[0][1] = 0.0;
    wgt[1][1] = 1.0;
    wgt[2][1] = 0.0;
    wgt[0][2] = 0.5;
    wgt[1][2] = 0.5;
    wgt[2][2] = 0.5;

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 5) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;

      mag = condition_number(p0,p1,p3,p4,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      mag = condition_number(p1,p2,p0,p4,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      mag = condition_number(p2,p3,p1,p4,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      mag = condition_number(p3,p0,p2,p4,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid corner condition number = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage pyramid corner condition number = %lg",favg);
      fprintf(out_f,"\nMaximum pyramid corner condition number = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnpri > 0)
  {
    wgt[0][0] = 1.0;
    wgt[1][0] = 0.0;
    wgt[2][0] = 0.0;
    wgt[0][1] = 0.5;
    wgt[1][1] = 0.5*sqrt(3);
    wgt[2][1] = 0.0;
    wgt[0][2] = 0.0;
    wgt[1][2] = 0.0;
    wgt[2][2] = 1.0;

    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 6) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;

      mag = condition_number(p0,p1,p2,p3,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      mag = condition_number(p1,p2,p0,p4,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      mag = condition_number(p2,p0,p1,p5,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      mag = condition_number(p3,p5,p4,p0,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
      mag = condition_number(p4,p3,p5,p1,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p4;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p4;
      }
      mag = condition_number(p5,p4,p3,p2,wgt);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p5;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p5;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism corner condition number = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage prism corner condition number = %lg",favg);
      fprintf(out_f,"\nMaximum prism corner condition number = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gnhex > 0)
  {
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 8) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      n6 = element[e].node_list->list[6];
      n7 = element[e].node_list->list[7];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;
      p6 = node[n6].vert;
      p7 = node[n7].vert;

      mag = tetrahedral_condition_number(p0,p1,p3,p4,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      mag = tetrahedral_condition_number(p1,p2,p0,p5,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      mag = tetrahedral_condition_number(p2,p3,p1,p6,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      mag = tetrahedral_condition_number(p3,p0,p2,p7,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
      mag = tetrahedral_condition_number(p4,p7,p5,p0,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p4;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p4;
      }
      mag = tetrahedral_condition_number(p5,p4,p6,p1,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p5;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p5;
      }
      mag = tetrahedral_condition_number(p6,p5,p7,p2,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p6;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p6;
      }
      mag = tetrahedral_condition_number(p7,p6,p4,p3,0);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p7;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p7;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral corner condition number = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage hexahedral corner condition number = %lg",favg);
      fprintf(out_f,"\nMaximum hexahedral corner condition number = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
    }
  }

  if (gntet > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 4) continue;
      if (element[e].me.proc != my_rank) continue;

      count++;
      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      cg = (p0+p1+p2+p3)/4.0;
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p2);
      v3 = Vector(p0,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      if (mag < 0.0)
        neg++;
      else
        pos++;

      if (mag < fmin)
      {
        fmin = mag;
        pmin = cg;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = cg;
      }
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
      i = neg;
      MPI_Allreduce(&i,&neg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = neg_skw;
      MPI_Allreduce(&i,&neg_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos_skw;
      MPI_Allreduce(&i,&pos_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos;
      MPI_Allreduce(&i,&pos,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral Jacobian = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage tetrahedral Jacobian = %lg",favg);
      fprintf(out_f,"\nMaximum tetrahedral Jacobian = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
      fprintf(out_f,"\nTetrahedral Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
    }
  }

  if (gnpyr > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 5) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;

      jmin = 1.0e20;
      javg = 0.0;
      jmax = -1.0e20;
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p3);
      v3 = Vector(p0,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
      javg /= 4.0;
      if (jmax < 0.0)
        neg++;
      else if (javg < 0.0 && jmax > 0.0)
        neg_skw++;
      else if (javg > 0.0 && jmin < 0.0)
        pos_skw++;
      else
        pos++;
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
      i = neg;
      MPI_Allreduce(&i,&neg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = neg_skw;
      MPI_Allreduce(&i,&neg_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos_skw;
      MPI_Allreduce(&i,&pos_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos;
      MPI_Allreduce(&i,&pos,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid Jacobian = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage pyramid Jacobian = %lg",favg);
      fprintf(out_f,"\nMaximum pyramid Jacobian = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
      fprintf(out_f,"\nPyramid Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
    }
  }

  if (gnpri > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 6) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;

      jmin = 1.0e20;
      javg = 0.0;
      jmax = -1.0e20;
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p2);
      v3 = Vector(p0,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      v1 = Vector(p2,p0);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      v1 = Vector(p3,p5);
      v2 = Vector(p3,p4);
      v3 = Vector(p3,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
      v1 = Vector(p4,p3);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p4;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p4;
      }
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p3);
      v3 = Vector(p5,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p5;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p5;
      }
      javg /= 6.0;
      if (jmax < 0.0)
        neg++;
      else if (javg < 0.0 && jmax > 0.0)
        neg_skw++;
      else if (javg > 0.0 && jmin < 0.0)
        pos_skw++;
      else
        pos++;
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
      i = neg;
      MPI_Allreduce(&i,&neg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = neg_skw;
      MPI_Allreduce(&i,&neg_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos_skw;
      MPI_Allreduce(&i,&pos_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos;
      MPI_Allreduce(&i,&pos,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism Jacobian = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage prism Jacobian = %lg",favg);
      fprintf(out_f,"\nMaximum prism Jacobian = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
      fprintf(out_f,"\nPrism Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
    }
  }

  if (gnhex > 0)
  {
    neg = neg_skw = pos_skw = pos = 0;
    fmin = 1.0e20;
    favg = 0.0;
    fmax = -1.0e20;
    count = 0;
    for (e=0; e < nelems; e++)
    {
      if (element_type(e) != 8) continue;
      if (element[e].me.proc != my_rank) continue;

      n0 = element[e].node_list->list[0];
      n1 = element[e].node_list->list[1];
      n2 = element[e].node_list->list[2];
      n3 = element[e].node_list->list[3];
      n4 = element[e].node_list->list[4];
      n5 = element[e].node_list->list[5];
      n6 = element[e].node_list->list[6];
      n7 = element[e].node_list->list[7];
      p0 = node[n0].vert;
      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;
      p4 = node[n4].vert;
      p5 = node[n5].vert;
      p6 = node[n6].vert;
      p7 = node[n7].vert;

      jmin = 1.0e20;
      javg = 0.0;
      jmax = -1.0e20;
      v1 = Vector(p0,p1);
      v2 = Vector(p0,p3);
      v3 = Vector(p0,p4);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p0;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p0;
      }
      v1 = Vector(p1,p2);
      v2 = Vector(p1,p0);
      v3 = Vector(p1,p5);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p1;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p1;
      }
      v1 = Vector(p2,p3);
      v2 = Vector(p2,p1);
      v3 = Vector(p2,p6);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p2;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p2;
      }
      v1 = Vector(p3,p0);
      v2 = Vector(p3,p2);
      v3 = Vector(p3,p7);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p3;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p3;
      }
      v1 = Vector(p4,p7);
      v2 = Vector(p4,p5);
      v3 = Vector(p4,p0);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p4;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p4;
      }
      v1 = Vector(p5,p4);
      v2 = Vector(p5,p6);
      v3 = Vector(p5,p1);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p5;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p5;
      }
      v1 = Vector(p6,p5);
      v2 = Vector(p6,p7);
      v3 = Vector(p6,p2);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p6;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p6;
      }
      v1 = Vector(p7,p6);
      v2 = Vector(p7,p4);
      v3 = Vector(p7,p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();
      mag = scalar_triple_product(v1,v2,v3);
      javg += mag;
      jmin = MIN(jmin,mag);
      jmax = MAX(jmax,mag);
      count++;
      if (mag < fmin)
      {
        fmin = mag;
        pmin = p7;
      }
      favg += mag;
      if (mag > fmax)
      {
        fmax = mag;
        pmax = p7;
      }
      javg /= 8.0;
      if (jmax < 0.0)
        neg++;
      else if (javg < 0.0 && jmax > 0.0)
        neg_skw++;
      else if (javg > 0.0 && jmin < 0.0)
        pos_skw++;
      else
        pos++;
    }
    favg /= MAX(1,count);

    if (my_rank == 0) fprintf(out_f,"\n");

    if (num_procs > 1)
    {
#ifdef PARALLEL
      g_in.val = fmin;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmin,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmin[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmin,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmin[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }

      favg *= count;
      MPI_Allreduce(&favg,&globald,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
      MPI_Allreduce(&count,&globali,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      favg = globald/globali;

      g_in.val = fmax;
      g_in.proc = my_rank;
      MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        if (g_out.proc != my_rank)
        {
          MPI_Recv(&fmax,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&pmax[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        }
      } else
      {
        if (g_out.proc == my_rank)
        {
          MPI_Send(&fmax,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&pmax[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
      i = neg;
      MPI_Allreduce(&i,&neg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = neg_skw;
      MPI_Allreduce(&i,&neg_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos_skw;
      MPI_Allreduce(&i,&pos_skw,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
      i = pos;
      MPI_Allreduce(&i,&pos,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral Jacobian = %lg, @ (%g, %g, %g)",fmin,pmin[0],pmin[1],pmin[2]);
      fprintf(out_f,"\nAverage hexahedral Jacobian = %lg",favg);
      fprintf(out_f,"\nMaximum hexahedral Jacobian = %lg, @ (%g, %g, %g)",fmax,pmax[0],pmax[1],pmax[2]);
      fprintf(out_f,"\nHexahedral Jacobian Statistics:");
      fprintf(out_f,"\n Negative        = %d",neg);
      fprintf(out_f,"\n Negative skewed = %d",neg_skw);
      fprintf(out_f,"\n Positive skewed = %d",pos_skw);
      fprintf(out_f,"\n Positive        = %d",pos);
    }
  }
  if (my_rank == 0) fprintf(out_f,"\n");

  return;
}
