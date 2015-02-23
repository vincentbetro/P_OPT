#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "geometry.h"
#include "smooth.h"
#include "Util.h"
#include "List.h"
#include "Bnormal.h"
#include "svdcmp.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

double pi = 4.0*atan(1.0);

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

int Smesh_obj::create_compress_row_storage(int **ia, int **iau, int **ja)
{
  List **nhash;
  int i, j, n;
  int mdim=0;

  // create node-to-node hash list
  nhash = new List*[nn];
  for (n=0; n < nn; n++)
    nhash[n] = new List();

  create_node_to_node(nhash);

  // create compressed row storage pointers
  *ia = new int[nn+1];
  for (n=0; n < nn; n++)
  {
    nhash[n]->Check_List(n);
    (*ia)[n] = mdim;
    mdim += nhash[n]->max;
  }
  (*ia)[nn] = mdim;
  *iau = new int[nn];
  *ja = new int[mdim];
  for (mdim=n=0; n < nn; n++)
  {
    for (i=0; i < nhash[n]->max; i++)
    {
      j = nhash[n]->list[i];
      if (j == n)
        (*iau)[n]=mdim;
      (*ja)[mdim++] = j;
    }
  }

  for (n=0; n < nn; n++)
    delete nhash[n];
  delete[] nhash;

  return(mdim);
}

int Smesh_obj::create_adjacency(int **ia, int **ja)
{
  List **nhash;
  int i, j, n;
  int mdim=0;

  // create node-to-node hash list
  nhash = new List*[nn];
  for (n=0; n < nn; n++)
    nhash[n] = new List();

  create_node_to_node(nhash);

  // create compressed row storage pointers
  *ia = new int[nn+1];
  for (n=0; n < nn; n++)
  {
    //nhash[n]->Check_List(n);
    (*ia)[n] = mdim;
    mdim += nhash[n]->max;
  }
  (*ia)[nn] = mdim;
  *ja = new int[mdim];
  for (mdim=n=0; n < nn; n++)
  {
    for (i=0; i < nhash[n]->max; i++)
    {
      j = nhash[n]->list[i];
      (*ja)[mdim++] = j;
    }
  }

  for (n=0; n < nn; n++)
    delete nhash[n];
  delete[] nhash;

  return(mdim);
}

void Smesh_obj::create_node_to_node(List **nhash)
{
  int c, n, n0, n1, n2, n3, n4, n5, n6, n7;

  for (n=0; n < nn; n++)
    nhash[n]->Redimension(0);

  for (c=0; c < ntet; c++)
  {
    n0 = tet_n[c][0];
    n1 = tet_n[c][1];
    n2 = tet_n[c][2];
    n3 = tet_n[c][3];
    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0)
      {
      fprintf(stderr,"\nYou have an index of %d %d %d %d in tet %d\n",n0,n1,n2,n3,c);
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
      }
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n2);
    nhash[n0]->Check_List(n3);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n3);
    nhash[n2]->Check_List(n0);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n1);
    nhash[n3]->Check_List(n2);
  }

  for (c=0; c < npyr; c++)
  {
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];

    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
      {
      fprintf(stderr,"\nYou have an index of %d %d %d %d %d in pyr %d\n",n0,n1,n2,n3,n4,c);
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
      }

    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n4);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n4);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n2);
    nhash[n4]->Check_List(n3);
  }

  for (c=0; c < npri; c++)
  {
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];

    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0)
      {
      fprintf(stderr,"\nYou have an index of %d %d %d %d %d %d in pri %d\n",n0,n1,n2,n3,n4,n5,c);
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
      }

    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n2);
    nhash[n0]->Check_List(n3);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n0);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n5);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n4);
    nhash[n3]->Check_List(n5);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n3);
    nhash[n4]->Check_List(n5);
    nhash[n5]->Check_List(n2);
    nhash[n5]->Check_List(n3);
    nhash[n5]->Check_List(n4);
  }

  for (c=0; c < nhex; c++)
  {
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];

    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0 || n6 < 0 || n7 < 0)
      {
      fprintf(stderr,"\nYou have an index of %d %d %d %d %d %d %d %d in hex %d\n",n0,n1,n2,n3,n4,n5,n6,n7,c);
      fflush(stderr);
#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,0);
#endif
      exit(0);
      }

    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n5);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n6);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n7);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n5);
    nhash[n4]->Check_List(n7);
    nhash[n5]->Check_List(n1);
    nhash[n5]->Check_List(n4);
    nhash[n5]->Check_List(n6);
    nhash[n6]->Check_List(n2);
    nhash[n6]->Check_List(n5);
    nhash[n6]->Check_List(n7);
    nhash[n7]->Check_List(n3);
    nhash[n7]->Check_List(n4);
    nhash[n7]->Check_List(n6);
  }

}

int Smesh_obj::critical_points(geometry *geom, int **cpn)
{
  int *tag;
  int g, i, j, k, n, ncp, n0, n1, n2, n3;

  tag = new int[nn];
  
  // identify critical points in mesh
  ncp = 0;
  if (nb == geom->ngb)
  {
    ncp = geom->ncp;
    *cpn = new int[ncp];
    for (k=0; k < ncp; k++)
    {
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      (*cpn)[k] = -1;

      for (n=0; n < nn; n++)
        tag[n] = 0;
      for (j=0; j < geom->c_point[k].nmat; j++)
      {
        g = geom->c_point[k].mat[j]-1;
        for (i=0; i < nt[g]; i++)
        {
          n0 = t_n[g][i][0];
          n1 = t_n[g][i][1];
          n2 = t_n[g][i][2];
          if (tag[n0] == j)
            tag[n0] = j+1;
          if (tag[n1] == j)
            tag[n1] = j+1;
          if (tag[n2] == j)
            tag[n2] = j+1;
        }
        for (i=0; i < nq[g]; i++)
        {
          n0 = q_n[g][i][0];
          n1 = q_n[g][i][1];
          n2 = q_n[g][i][2];
          n3 = q_n[g][i][3];
          if (tag[n0] == j)
            tag[n0] = j+1;
          if (tag[n1] == j)
            tag[n1] = j+1;
          if (tag[n2] == j)
            tag[n2] = j+1;
          if (tag[n3] == j)
            tag[n3] = j+1;
        }
      }

      double dsmn = 1.0e20;
      i = -1;
      for (n=0; n < nn; n++)
      {
#ifdef PARALLEL
        if (pmap[n][1] == my_rank && tag[n] == geom->c_point[k].nmat)
#else
        if (tag[n] == geom->c_point[k].nmat)
#endif
        {
          double ds = distance(geom->g_vert[geom->c_point[k].index],nodep[n]);
          if (ds < dsmn)
          {
            dsmn = ds;
            i = n;
          }
        }
      }

#ifdef PARALLEL
      MPI_Status single_status;
      //use struct to get MPI_MINLOC
      struct 
      { 
        double val; 
        int proc; 
      } dsmn_g_in, dsmn_g_out; 

      dsmn_g_in.val = dsmn;
      dsmn_g_in.proc = my_rank;
      dsmn_g_out.val = dsmn;
      dsmn_g_out.proc = my_rank;

      MPI_Allreduce(&dsmn_g_in,&dsmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);

      Point gp, np;
      if (my_rank == 0)
      {
        if (dsmn_g_out.proc != my_rank)
        {
          MPI_Recv(&i,1,MPI_INT,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&gp[0],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&gp[1],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&gp[2],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&np[0],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&np[1],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
          MPI_Recv(&np[2],1,MPI_DOUBLE,dsmn_g_out.proc,dsmn_g_out.proc,MPI_COMM_WORLD, &single_status);
        } else
        {
          (*cpn)[k] = i;
          i = pmap[i][0];
          gp = geom->g_vert[geom->c_point[k].index];
          np = nodep[(*cpn)[k]];
        }
        fprintf(out_f,"\nCritical point %d node index = %d",k+1,i);
        fprintf(out_f,"\n  Geometry location = ( %lg, %lg, %lg )",gp[0],gp[1],gp[2]);
        fprintf(out_f,"\n      Mesh location = ( %lg, %lg, %lg )",np[0],np[1],np[2]);
        fflush(out_f);
      } else
      {
        if (dsmn_g_out.proc == my_rank)
        {
          (*cpn)[k] = i;
          i = pmap[i][0];
          gp = geom->g_vert[geom->c_point[k].index];
          np = nodep[(*cpn)[k]];
          MPI_Send(&i,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&gp[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&gp[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&gp[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&np[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&np[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
          MPI_Send(&np[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        }
      }
#else
      (*cpn)[k] = i;
      fprintf(out_f,"\nCritical point %d node index = %d",k+1,(*cpn)[k]);
      fprintf(out_f,"\n  Geometry location = ( %lg, %lg, %lg )",
              geom->g_vert[geom->c_point[k].index][0],
              geom->g_vert[geom->c_point[k].index][1],
              geom->g_vert[geom->c_point[k].index][2]);
      fprintf(out_f,"\n      Mesh location = ( %lg, %lg, %lg )",
              nodep[(*cpn)[k]][0],nodep[(*cpn)[k]][1],nodep[(*cpn)[k]][2]);
      fflush(out_f);
#endif
    }
  }

  delete[] tag;

  return(ncp);
}

//no need to parallelie, since just checking proc by proc for problems, and even checks ghost cells

int Smesh_obj::check_volumes()
{
  int c, i, k, n, n0, n1, n2, n3, n4, n5, n6, n7, success, source, dest, local, tag;
  double vol, vmin, gdouble;
  Point cp, p0, p1, p2, p3, p4, p5, p6, p7;
  int gntet, gnpyr, gnpri, gnhex;

#ifdef PARALLEL
  MPI_Status status;
#endif

  //use struct to get MPI_MINLOC & MPI_MAXLOC
  struct 
  { 
    double val; 
    int proc; 
  } cmn_g_in, cmn_g_out, cmx_g_in, cmx_g_out; 

  success = 1;

  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;
#ifdef PARALLEL
  if (num_procs > 1)
  {
    MPI_Barrier(MPI_COMM_WORLD);

    gntet=local=0;
    for (n=0; n < ntet; n++)
      local = MAX(local,tet_map[n]+1);
    gntet = local;
    MPI_Allreduce(&local,&gntet,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    gnpyr=local=0;
    for (n=0; n < npyr; n++)
      local = MAX(local,pyr_map[n]+1);
    gnpyr = local;
    MPI_Allreduce(&local,&gnpyr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    gnpri=local=0;
    for (n=0; n < npri; n++)
      local = MAX(local,pri_map[n]+1);
    gnpri = local;
    MPI_Allreduce(&local,&gnpri,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

    gnhex=local=0;
    for (n=0; n < nhex; n++)
      local = MAX(local,hex_map[n]+1);
    gnhex = local;
    MPI_Allreduce(&local,&gnhex,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  }
#endif

  // check cell volumes

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int cmin;
  vmin = 1.0e-40;
  if (gntet > 0)
  {
    for (k=c=0; c < ntet; c++)
    {
      n0 = tet_n[c][0];
      n1 = tet_n[c][1];
      n2 = tet_n[c][2];
      n3 = tet_n[c][3];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      vol = tetrahedral_volume(p0,p1,p2,p3);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          cmin = c;
          vmin = vol;
          cp = (p0+p1+p2+p3)/4.0;
        }
        success=0;
        k++;
      }
    }
    i = k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum tetrahedral vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      if (num_procs == 1)
      {
        fprintf(out_f,"\n  tet = %d",cmin);
        c = cmin;
        n0 = tet_n[c][0];
        n1 = tet_n[c][1];
        n2 = tet_n[c][2];
        n3 = tet_n[c][3];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        fprintf(out_f,"\n Node 0 = %d",n0); p0.print(out_f);
        fprintf(out_f,"\n Node 1 = %d",n1); p1.print(out_f);
        fprintf(out_f,"\n Node 2 = %d",n2); p2.print(out_f);
        fprintf(out_f,"\n Node 3 = %d",n3); p3.print(out_f);
      }
      fprintf(out_f,"\nTotal number of negative tetrahedral volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpyr > 0)
  {
    for (k=c=0; c < npyr; c++)
    {
      n0 = pyr_n[c][0];
      n1 = pyr_n[c][1];
      n2 = pyr_n[c][2];
      n3 = pyr_n[c][3];
      n4 = pyr_n[c][4];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      vol = pyramid_volume(p0,p1,p2,p3,p4);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          cmin = c;
          vmin = vol;
          cp = (p0+p1+p2+p3+p4)/5.0;
        }
        success=0;
        k++;
      }
    }
    i = k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum pyramid vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      if (num_procs == 1)
      {
        fprintf(out_f,"\n  pyr = %d",cmin);
        c = cmin;
        n0 = pyr_n[c][0];
        n1 = pyr_n[c][1];
        n2 = pyr_n[c][2];
        n3 = pyr_n[c][3];
        n4 = pyr_n[c][4];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        fprintf(out_f,"\n Node 0 = %d",n0); p0.print(out_f);
        fprintf(out_f,"\n Node 1 = %d",n1); p1.print(out_f);
        fprintf(out_f,"\n Node 2 = %d",n2); p2.print(out_f);
        fprintf(out_f,"\n Node 3 = %d",n3); p3.print(out_f);
        fprintf(out_f,"\n Node 4 = %d",n4); p4.print(out_f);
      }
      fprintf(out_f,"\nTotal number of negative pyramid volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnpri > 0)
  {
    for (k=c=0; c < npri; c++)
    {
      n0 = pri_n[c][0];
      n1 = pri_n[c][1];
      n2 = pri_n[c][2];
      n3 = pri_n[c][3];
      n4 = pri_n[c][4];
      n5 = pri_n[c][5];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      vol = prism_volume(p0,p1,p2,p3,p4,p5);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          cmin = c;
          vmin = vol;
          cp = (p0+p1+p2+p3+p4+p5)/6.0;
        }
        success=0;
        k++;
      }
    }
    i = k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum prism vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      if (num_procs == 1)
      {
        fprintf(out_f,"\n  pri = %d",cmin);
        c = cmin;
        n0 = pri_n[c][0];
        n1 = pri_n[c][1];
        n2 = pri_n[c][2];
        n3 = pri_n[c][3];
        n4 = pri_n[c][4];
        n5 = pri_n[c][5];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        p5 = nodep[n5];
        fprintf(out_f,"\n Node 0 = %d",n0); p0.print(out_f);
        fprintf(out_f,"\n Node 1 = %d",n1); p1.print(out_f);
        fprintf(out_f,"\n Node 2 = %d",n2); p2.print(out_f);
        fprintf(out_f,"\n Node 3 = %d",n3); p3.print(out_f);
        fprintf(out_f,"\n Node 4 = %d",n4); p4.print(out_f);
        fprintf(out_f,"\n Node 5 = %d",n5); p5.print(out_f);
      }
      fprintf(out_f,"\nTotal number of negative prism volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  vmin = 1.0e-40;
  if (gnhex > 0)
  {
    for (k=c=0; c < nhex; c++)
    {
      n0 = hex_n[c][0];
      n1 = hex_n[c][1];
      n2 = hex_n[c][2];
      n3 = hex_n[c][3];
      n4 = hex_n[c][4];
      n5 = hex_n[c][5];
      n6 = hex_n[c][6];
      n7 = hex_n[c][7];
      p0 = nodep[n0];
      p1 = nodep[n1];
      p2 = nodep[n2];
      p3 = nodep[n3];
      p4 = nodep[n4];
      p5 = nodep[n5];
      p6 = nodep[n6];
      p7 = nodep[n7];
      vol = hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7);
      if (vol < 1.0e-40)
      {
        if (vol < vmin)
        {
          cmin = c;
          vmin = vol;
          cp = (p0+p1+p2+p3+p4+p5+p6+p7)/8.0;
        }
        success=0;
        k++;
      }
    }
    i = k;
#ifdef PARALLEL
    MPI_Allreduce(&k,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&vmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    cmn_g_in.val = vmin;
    cmn_g_in.proc = my_rank;
    vmin = gdouble;
    MPI_Allreduce(&cmn_g_in,&cmn_g_out,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    if (cmn_g_out.proc != 0)
    {
      if (my_rank == cmn_g_out.proc)
      {
        tag = my_rank;
        dest = 0;
        MPI_Send(&cp[0], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(&cp[2], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        tag = source = cmn_g_out.proc;
        MPI_Recv(&cp[0], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&cp[2], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      }
    }
#endif
    if (i > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nMinimum hexahedral vol = %lg, x,y,z = %g, %g, %g\n",vmin,cp[0],cp[1],cp[2]);
      if (num_procs == 1)
      {
        fprintf(out_f,"\n  hex = %d",cmin);
        c = cmin;
        n0 = hex_n[c][0];
        n1 = hex_n[c][1];
        n2 = hex_n[c][2];
        n3 = hex_n[c][3];
        n4 = hex_n[c][4];
        n5 = hex_n[c][5];
        n6 = hex_n[c][6];
        n7 = hex_n[c][7];
        p0 = nodep[n0];
        p1 = nodep[n1];
        p2 = nodep[n2];
        p3 = nodep[n3];
        p4 = nodep[n4];
        p5 = nodep[n5];
        p6 = nodep[n6];
        p7 = nodep[n7];
        fprintf(out_f,"\n Node 0 = %d",n0); p0.print(out_f);
        fprintf(out_f,"\n Node 1 = %d",n1); p1.print(out_f);
        fprintf(out_f,"\n Node 2 = %d",n2); p2.print(out_f);
        fprintf(out_f,"\n Node 3 = %d",n3); p3.print(out_f);
        fprintf(out_f,"\n Node 4 = %d",n4); p4.print(out_f);
        fprintf(out_f,"\n Node 5 = %d",n5); p5.print(out_f);
        fprintf(out_f,"\n Node 6 = %d",n6); p6.print(out_f);
        fprintf(out_f,"\n Node 7 = %d",n7); p7.print(out_f);
      }
      fprintf(out_f,"\nTotal number of negative hexahedral volumes = %d",i);
      fflush(out_f);
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(&success,&i,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
#endif
  success = i;
  return(success);
}

