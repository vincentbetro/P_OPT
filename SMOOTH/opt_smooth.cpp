#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "geometry.h"
#include "Point.h"
#include "Vector.h"
#include "smooth.h"
#include "Util.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/
extern int rotate;
extern int blend_exponent;
 
double pi = 4.0*atan(1.0);
double Jzero = 0.0;
bool Exp_cost = true;

int POLYMESH::critical_points(geometry *geom, int **cpn)
{
  int *tag;
  int f, g, i, j, k, m, n, ncp;

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
        for (i=0; i < boundary[g].polygon_list->max; i++)
        {
          f = abs(boundary[g].polygon_list->list[i]);
          for (m=0; m < face[f].node_list->max; m++)
          {
            n = face[f].node_list->list[m];
            if (tag[n] == j)
              tag[n] = j+1;
          }
        }
      }

      double dsmn = 1.0e20;
      i = -1;
      for (n=0; n < nn; n++)
      {
#ifdef PARALLEL
        if (node[k].me.proc == my_rank && tag[n] == geom->c_point[k].nmat)
#else
        if (tag[n] == geom->c_point[k].nmat)
#endif
        {
          double ds = distance(geom->g_vert[geom->c_point[k].index],node[n].vert);
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
          i = node[i].me.global;
          gp = geom->g_vert[geom->c_point[k].index];
          np = node[(*cpn)[k]].vert;
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
          i = node[i].me.global;
          gp = geom->g_vert[geom->c_point[k].index];
          np = node[(*cpn)[k]].vert;
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
              node[(*cpn)[k]].vert[0],node[(*cpn)[k]].vert[1],node[(*cpn)[k]].vert[2]);
      fflush(out_f);
#endif
    }
  }

  delete[] tag;

  return(ncp);
}

int POLYMESH::opt_smooth(geometry *geom, int nsmoo, int mvbnd, int mvint, int sflag, int eflag, double threshold)
{
  int b, f, g, i, j, k, l, m, n, ncp;
  int sweep, nmx;
  int *tag;
  int *cpn;
  List **blist, *bnode, **flist;
  double ds, dsmx, mag;
  Point pmx, del;
  Point p0, p1, p2, p3, p4;
  Vector v1, v2, v3, norm;
  bool surf_cost = true;

  int success = 1;

  if (mvbnd)
    if (my_rank == 0) fprintf(out_f,"\nMoving boundary nodes!\n");
  if (mvint < 0)
    if (my_rank == 0) fprintf(out_f,"\nMoving all but %d interior node layers!\n",-mvint);
  if (mvint > 0)
    if (my_rank == 0) fprintf(out_f,"\nMoving %d interior node layers!\n",mvint);
  fflush(out_f);

  // identify critical points in mesh
  ncp = critical_points(geom,&cpn);

  // create face list for boundary nodes
  flist = new List*[nn];
  for (n=0; n < nn; n++)
    flist[n] = new List();

  bnode = new List();

  int *bindex = new int[nn];
  for (n=0; n < nn; n++)
    bindex[n] = -1;

  for (b=0; b < nb; b++)
  {
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = boundary[b].polygon_list->list[i];
      for (j=0; j < face[abs(f)].node_list->max; j++)
      {
        n = face[abs(f)].node_list->list[j];
        bindex[n] = 1;
        flist[n]->Check_List(f);
      }
    }
  }
  for (n=0; n < nn; n++)
    if (bindex[n] > 0)
    {
      bindex[n] = bnode->max; // index into bnode list
      bnode->Add_To_List(n);
    }

  blist = new List*[bnode->max];
  for (n=0; n < bnode->max; n++)
    blist[n] = new List();

  for (b=0; b < nb; b++)
  {
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = abs(boundary[b].polygon_list->list[i]);
      for (j=0; j < face[f].node_list->max; j++)
      {
        n = face[f].node_list->list[j];
        blist[bindex[n]]->Check_List(b+1);
      }
    }
  }

  // set tag to 0 for interior
  // set to -1 for boundary nodes not on edge or critical point
  // set to -2 for boundary nodes on edge
  // set to -3 for boundary nodes at critical point

  tag = new int[nn];
  for (n=0; n < nn; n++)
  {
    i = bindex[n];
    if (i < 0) tag[n] = 0;
    if (i >= 0 && blist[i]->max == 1) tag[n] = -1;
    if (i >= 0 && blist[i]->max == 2) tag[n] = -2;
    if (i >= 0 && blist[i]->max > 2) tag[n] = -3;
  }

  for (k=n=0; n < nn; n++)
    if (bindex[n] >= 0)
      k++;
#ifdef PARALLEL
    i = k;
    MPI_Allreduce(&i,&k,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
  if (my_rank == 0) fprintf(out_f,"\nNumber of boundary nodes = %d",k);
  fflush(out_f);


  // set tag for interior nodes by marching away from walls
  int globalg;
  globalg = g = 1;
  m = 0;
  while (globalg > 0)
  {
    exchange_node_int(tag);
    g = 0;
    m++;
    for (n=0; n < nn; n++)
    {
      if (node[n].me.proc != my_rank || tag[n] != 0)
        continue;
      for (i=0; i < node[n].nabor->max && tag[n] == 0; i++)
      {
        j = node[n].nabor->list[i];
        if (tag[j] != 0 && tag[j] < m)
        {
          tag[n] = m;
          g++;
        }
      }
    }
    globalg = g;
#ifdef PARALLEL
    MPI_Allreduce(&g,&globalg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);        
#endif
    if (my_rank == 0) fprintf(out_f,"\nNumber of nodes tagged with %d = %d",m,globalg);
    fflush(out_f);
  }
  if (my_rank == 0) fprintf(out_f,"\nMaximum node tag number = %d",m-1);
  fflush(out_f);

  // project boundary nodes to surface
  ds = 0.0;
  for (l=0; l < bnode->max; l++)
  {
    n = bnode->list[l];
    if (node[n].me.proc != my_rank)
      continue;
    bool bflag = false;
    if (blist[l]->max < 3)
    {
      bflag = true;
      for (j=0; j < blist[l]->max && bflag; j++)
        if (abs(geom->layers[blist[l]->list[j]]) != 1)
          bflag = false;
    }
    if (bflag)
    {
      p1 = node[n].vert;
      if (blist[l]->max == 1)
      {
        norm = Vector(0.0,0.0,0.0);
        p2 = geom->closest(p1,blist[l]->list[0],norm);
        node[n].vert = p2;
        ds = MAX(ds,distance(p1,p2));
      }
      if (blist[l]->max == 2)
      {
        p2 = geom->closest_chain(p1,blist[l]->list[0],blist[l]->list[1]);
        node[n].vert = p2;
        ds = MAX(ds,distance(p1,p2));
      }
    }
  }
  dsmx = ds;
#ifdef PARALLEL
  MPI_Allreduce(&ds,&dsmx,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  if (dsmx > 1.0e-12 && my_rank == 0)
  {
    fprintf(out_f,"\nMaximum boundary displacement = %lg",dsmx);
    fflush(out_f);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Check for movement restriction flags
  int movex=1, movey=1, movez=1;

  FILE *fp;
  if ((fp = fopen("P_OPT.move_restrict","r")) != NULL)
  {
    fscanf(fp,"%d %d %d",&movex,&movey,&movez);
    if (my_rank == 0) fprintf(out_f,"\nInput values for movex, movey, movez = %d, %d, %d",movex,movey,movez);
    fclose(fp);
  }

  Point box[2];

  box[0] = geom->min_point();
  box[1] = geom->max_point();

  box[0] -= Point(1.0,1.0,1.0);
  box[1] += Point(1.0,1.0,1.0);

  if (my_rank == 0)
  {
    fprintf(out_f,"\nDefault values for bounding box:");
    fprintf(out_f,"\n   low corner = (%lf, %lf, %lf)",box[0][0],box[0][1],box[0][2]);
    fprintf(out_f,"\n  high corner = (%lf, %lf, %lf)",box[1][0],box[1][1],box[1][2]);
    fflush(out_f);
  }

  // Check for movement restriction box
  if ((fp = fopen("P_OPT.move_restrict_box","r")) != NULL)
  {
    fscanf(fp,"%lf %lf %lf",&(box[0][0]),&(box[0][1]),&(box[0][2]));
    fscanf(fp,"%lf %lf %lf",&(box[1][0]),&(box[1][1]),&(box[1][2]));
    fclose(fp);

    Point lo, hi;
    for (i=0; i < 3; i++)
    {
      lo[i] = MIN(box[0][i],box[1][i]);
      hi[i] = MAX(box[0][i],box[1][i]);
    }
    box[0] = lo;
    box[1] = hi;

    if (my_rank == 0)
    {
      fprintf(out_f,"\nInput values for bounding box:");
      fprintf(out_f,"\n   low corner = (%lf, %lf, %lf)",box[0][0],box[0][1],box[0][2]);
      fprintf(out_f,"\n  high corner = (%lf, %lf, %lf)",box[1][0],box[1][1],box[1][2]);
      fflush(out_f);
    }
  }

  // compute epsilon distance;
  double eps = 1.0e20;
  int kk;
  for (kk=n=0; n < nn; n++)
  {
    double dsmn = 0.0;
    k=0;
    for (i=0; i < node[n].nabor->max; i++)
    {
      j = node[n].nabor->list[i];
      ds = distance(node[n].vert,node[j].vert);
      dsmn += ds;
      if (ds < 1.0e-15) k=1;
    }
    if (k)
    {
      kk++;
      p0=Point(0.0,0.0,0.0);
      for (i=0; i < node[n].nabor->max; i++)
      {
        j = node[n].nabor->list[i];
        p0 += node[j].vert;
      }
      node[n].vert = p0/node[n].nabor->max;
    } else
      eps = MIN(eps,dsmn/node[n].nabor->max);
  }
#ifdef PARALLEL
  k = kk;
  MPI_Allreduce(&kk,&k,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (k > 0 && my_rank == 0)
    fprintf(out_f,"\n\nNumber of points repositioned = %d\n",k);

  eps = MAX(1.0e-15,eps*1.0e-03);
  mag = eps;
#ifdef PARALLEL
  MPI_Allreduce(&mag,&eps,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  if (my_rank == 0) fprintf(out_f,"\nepsilon = %g",eps);
  if (my_rank == 0) fprintf(out_f,"\ncost threshold = %g",threshold);

  // begin main smoothing loop
  Vector *pert = new Vector[nn];

  double cost = 0.0, cnew = 0.0;

  if (my_rank == 0) fprintf(out_f,"\n\n Iter   # exceed  # inverted  # perturb  Max Delta     Min Cost     Avg Cost     Max Cost    Max Cost X   Max Cost Y   Max Cost Z");
  if (my_rank == 0) fprintf(out_f,  "\n------ ---------- ---------- ---------- ------------ ------------ ------------ ------------ ------------ ------------ ------------");
  fflush(out_f);

  for (sweep=1; sweep <= MAX(1,abs(nsmoo)); sweep++)
  {

    for (n=0; n < nn; n++)
      pert[n] = Vector(0.0,0.0,0.0);

    dsmx = 0.0;
    double cavg= 0.0;
    double cmn = 1.0e20;
    double cmx = 0.0;
    nmx = -1;
    int knt = 0;
    int npert = 0;
    int exceed = 0;
    int neg = 0;
        
    int start, finish, inc;
    if (sweep % 2 == 0)
    {
      start = nn-1;
      finish = 0;
      inc = -1;
    } else
    {
      start = 0;
      finish = nn-1;
      inc = 1;
    }

    for (n=start; n != finish; n += inc)
    {
      if (node[n].me.proc != my_rank) continue;

      if (tag[n] == -3) continue;

      // check movement restriction within bounding box
      if (node[n].vert[0] < box[0][0] || node[n].vert[0] > box[1][0] ||
          node[n].vert[1] < box[0][1] || node[n].vert[1] > box[1][1] ||
          node[n].vert[2] < box[0][2] || node[n].vert[2] > box[1][2])
          continue;
  
      // compute perturbation amount for current node
      double pmag=1.0e20;
      for (i=0; i < node[n].nabor->max; i++)
      {
        j = node[n].nabor->list[i];
        pmag = MIN(pmag,distance(node[j].vert,node[n].vert));
      }
      pmag /= 20.0;
      pmag = MAX(pmag,eps);

      if (eflag == 0)
        Exp_cost = false;
      else
        Exp_cost = true;

      if (sflag == 1)
        surf_cost = true;
      else
        surf_cost = false;

      if (bindex[n] >= 0)
      {
        // boundary point
        if (mvbnd == 0) continue;

        // check for floating boundary
        bool bflag = true;
        if ((i=bindex[n]) >= 0)
        {
          if (blist[i]->max > 2) bflag = false;
          for (j=0; j < blist[i]->max && bflag; j++)
            if (geom->layers[blist[i]->list[j]] >= 0)
              bflag = false;
        }
        if (!bflag) continue;

      } else
      {
        // interior point
        if (mvint == 0) continue;
        if ((mvint > 0 && tag[n] > mvint) || (mvint < 0 && tag[n] < abs(mvint))) continue;
      }

      if (tag[n] < 0 && surf_cost)
        cost = cost_function(geom,n,flist[n],bindex,blist);
      else
        cost = cost_function(n);

      if (cost >= 1.0)
      {
        neg++;
        if (tag[n] < 0)
        {
          surf_cost = true;
          cost = cost_function(geom,n,flist[n],bindex,blist);
        }
      }

      cmn = MIN(cmn,cost);
      cavg += cost;
      if (cost > cmx)
      {
        cmx = cost;
        pmx = node[n].vert;
        nmx = n;
      }
      knt++;

      if (cost < threshold) continue;

      double pert_cost = cost;
      exceed++;

      int cube_dim = 1;
      cube_dim  = (cost >= 1.0) ? 3 : 1 ;
      for (k=-cube_dim; k <= cube_dim; k++)
      {
        for (j=-cube_dim; j <= cube_dim; j++)
        {
          for (i=-cube_dim; i <= cube_dim; i++)
          {
            if (i==0 && j==0 && k==0) continue;
            if (abs(i) < cube_dim && abs(j) < cube_dim && abs(k) < cube_dim) continue;

            Vector vec1 = Vector((double)i,(double)j,(double)k);
            vec1.normalize();
            vec1 *= pmag;
            double dx = vec1[0]*movex;
            double dy = vec1[1]*movey;
            double dz = vec1[2]*movez;
            Point original = node[n].vert;
            node[n].vert[0] += dx;
            node[n].vert[1] += dy;
            node[n].vert[2] += dz;
            if (tag[n] == -2)
            {
              l = bindex[n];
              p2 = geom->closest_chain(node[n].vert,blist[l]->list[0],blist[l]->list[1]);
              node[n].vert = p2;
              dx = node[n].vert[0]-original[0];
              dy = node[n].vert[1]-original[1];
              dz = node[n].vert[2]-original[2];
            }
            if (tag[n] == -1)
            {
              l = bindex[n];
              norm = Vector(0.0,0.0,0.0);
              p2 = geom->closest(node[n].vert,blist[l]->list[0],norm);
              node[n].vert = p2;
              dx = node[n].vert[0]-original[0];
              dy = node[n].vert[1]-original[1];
              dz = node[n].vert[2]-original[2];
            }

            if (tag[n] < 0 && surf_cost)
              cnew = cost_function(geom,n,flist[n],bindex,blist);
            else
              cnew = cost_function(n);

            if (cnew < pert_cost)
            {
              pert_cost = cnew;
              pert[n] = Vector(dx,dy,dz);
            }
            node[n].vert[0] -= dx;
            node[n].vert[1] -= dy;
            node[n].vert[2] -= dz;
          }
        }
      }

      del[0] = pert[n][0];
      del[1] = pert[n][1];
      del[2] = pert[n][2];
      ds=(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
      dsmx = MAX(dsmx,ds);
      if (pert[n].magnitude() > 1.0e-12) npert++;

      if (nsmoo > 0)
        node[n].vert += Point(pert[n][0],pert[n][1],pert[n][2]);
    }
    dsmx = sqrt(dsmx);

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status single_status;
    struct 
    { 
      double val; 
      int proc; 
    } g_in, g_out; 

    g_in.val = cmx;
    g_in.proc = my_rank;
    MPI_Allreduce(&g_in,&g_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    if (my_rank == 0)
    {
      if (g_out.proc != my_rank)
      {
        MPI_Recv(&cmx,1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&nmx,1,MPI_INT,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmx[0],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmx[1],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
        MPI_Recv(&pmx[2],1,MPI_DOUBLE,g_out.proc,g_out.proc,MPI_COMM_WORLD, &single_status);
      } else
      {
        nmx = node[nmx].me.global;
      }
    } else
    {
      if (g_out.proc == my_rank)
      {
        MPI_Send(&cmx,1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&node[nmx].me.global,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmx[0],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmx[1],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
        MPI_Send(&pmx[2],1,MPI_DOUBLE,0,my_rank,MPI_COMM_WORLD);
      }
    }

    mag = dsmx;
    MPI_Allreduce(&mag,&dsmx,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    mag = cmn;
    MPI_Allreduce(&mag,&cmn,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

    mag = cavg;
    MPI_Allreduce(&mag,&cavg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    i = exceed;
    MPI_Allreduce(&i,&exceed,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    i = knt;
    MPI_Allreduce(&i,&knt,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    i = neg;
    MPI_Allreduce(&i,&neg,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    i = npert;
    MPI_Allreduce(&i,&npert,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#endif

    cavg /= MAX(1,knt);

    if (my_rank == 0)
      fprintf(out_f,"\n%6d %10d %10d %10d %12.6e %12.6e %12.6e %12.6e %12.5e %12.5e %12.5e",
       sweep,exceed,neg,npert,dsmx,cmn,cavg,cmx,pmx[0],pmx[1],pmx[2]);

    fflush(out_f);

    if (num_procs > 1)
      exchange_node_vector(pert);

    if (nsmoo > 0)
    {
      for (n=0; n < nn; n++)
        if (node[n].me.proc != my_rank)
          node[n].vert += Point(pert[n][0],pert[n][1],pert[n][2]);
    }

    // Ensure critical points don't move
    if (mvbnd)
    {
      for (k=0; k < ncp; k++)
      {
        n = cpn[k];
        if (n < 0)
          continue;
        node[n].vert = geom->g_vert[geom->c_point[k].index];
      }
    }

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (npert == 0 && cmx < threshold)
      break;
  }

  //success = check_volumes();

  // clean up memory
  if (nb == geom->ngb && ncp > 0)
    delete[] cpn;
  if (bnode->max > 0)
  {
    for (n=0; n < bnode->max; n++)
      delete blist[n];
    delete[] blist;
    delete bnode;
  }
  for (n=0; n < nn; n++)
    delete flist[n];
  delete[] flist;
  delete[] bindex;
  delete[] pert;
  delete[] tag;

  return(success);
}

bool POLYMESH::node_in_element(int n, int e)
{
  bool flag = false;

  if (element[e].node_list->Is_In_List(n))
    flag = true;

  return(flag);
}

int POLYMESH::adjacent_nodes(int e, int m, List *nlist, List *flist)
{
  int f, j, k, n1, n2, n3;

  // create list of nodes in element e connected to m and wound properly
  // at the same time create a list of faces used to build node list
  nlist->max = 0;
  flist->max = 0;
  if (!node_in_element(m,e))
    return(0);

  bool done = false;
  bool changed = true;
  while (!done && changed)
  {
    changed = false;
    for (j=0; j < element[e].face_list->max && !done; j++)
    {
      f = element[e].face_list->list[j];
      if (!face[abs(f)].node_list->Is_In_List(m))
        continue;

      for (k=0; k < face[abs(f)].node_list->max; k++)
      {
        n1 = face[abs(f)].node_list->list[k];
        if (n1 != m)
          continue;
        if (k < face[abs(f)].node_list->max-1)
          n3 = face[abs(f)].node_list->list[k+1];
        else
          n3 = face[abs(f)].node_list->list[0];
        if (k == 0)
          n2 = face[abs(f)].node_list->list[face[abs(f)].node_list->max-1];
        else
          n2 = face[abs(f)].node_list->list[k-1];

        if (f < 0)
        {
          int tmp = n2;
          n2 = n3;
          n3 = tmp;
        }
        if (nlist->max == 0)
        {
          nlist->Add_To_List(n2);
          changed = true;
        }
        if (n2 == nlist->list[nlist->max-1])
        {
          if (n3 == nlist->list[0])
            done = true;
          else
            nlist->Add_To_List(n3);
          flist->Add_To_List(f);
          changed = true;
        }
      }
    }
  }
  if (!done)
    return(-1);
  else
    return(nlist->max);
}

double POLYMESH::cost_function(geometry *geom, int n, List *flist, int bindex[], List **blist)
{
  // compute cost based on surface elements
  double angle, cost, cmin, cmax, c, cn, J;
  int count, f, g, i, j, k, n1, n2, n3;
  int ii, jj, kk, i1, i2, i3, b1, b2, b3;
  Point cg, p1, p2, p3;
  Vector u, v, w;
  double wgt[2][2];

  cmax = cost = 0.0;
  cmin = 1.0e20;
  count = 0;

  for (i=0; i < flist->max; i++)
  {
    f = flist->list[i];

    angle = 180.0 - 360.0/face[abs(f)].node_list->max;
    angle = angle*pi/180.0;

    wgt[0][0]=1.0;
    wgt[1][0]=0.0;
    wgt[0][1]=cos(angle);
    wgt[1][1]=sin(angle);

    for (k=0; k < face[abs(f)].node_list->max; k++)
    {
      n1 = face[abs(f)].node_list->list[k];
      if (k < face[abs(f)].node_list->max-1)
        n2 = face[abs(f)].node_list->list[k+1];
      else
        n2 = face[abs(f)].node_list->list[0];
      if (k > 0)
        n3 = face[abs(f)].node_list->list[k-1];
      else
        n3 = face[abs(f)].node_list->list[face[abs(f)].node_list->max-1];

      if (n1 != n && n2 != n && n3 != n)
        continue;

      if (!Exp_cost && n1 != n)
        continue;

      p1 = node[n1].vert;
      p2 = node[n2].vert;
      p3 = node[n3].vert;

      cg = (p1+p2+p3)/3.0;
      Vector norm = Vector(0.0,0.0,0.0);
      i1 = bindex[n1];
      i2 = bindex[n2];
      i3 = bindex[n3];
      bool flag = false;
      for (ii=0; ii < blist[i1]->max && !flag; ii++)
      {
        b1 = blist[i1]->list[ii];
        for (jj=0; jj < blist[i2]->max && !flag; jj++)
        {
          b2 = blist[i2]->list[jj];
          if (b2 != b1) continue;
          for (kk=0; kk < blist[i3]->max && !flag; kk++)
          {
            b3 = blist[i3]->list[kk];
            if (b1 == b3)
              flag = true;
          }
        }
      }
      if (!flag)
        continue;
      
      cg = geom->closest(cg,b1,norm);

      u = Vector(p1,p2);
      v = Vector(p1,p3);
      if (f < 0)
        w = v % u;
      else
        w = u % v;
      
      J = w * norm;
      if (J < Jzero)
      {
        c = 1.0-J;
      } else
      {
        cn = condition_number(p1,p2,p3,wgt);
        c = (1.0-1.0/cn);
      }
      cost += c;
      count++;
      cmax = MAX(cmax,c);
    }
  }

  cost/=MAX(1,count);
  double ratio = pow(cmax,(double)blend_exponent);
  cost=ratio*cmax + MAX(0.0,(1.0-ratio))*cost;
  return (cost);
}

double POLYMESH::cost_function(int n)
{
  // compute cost based on volume elements
  double cost, cmin, cmax, c, J;
  int count, e, h, i, j, k, n1, n2, n3, f1, f2, f3;
  int jmax;
  Point p0, p1, p2, p3;
  Vector v1, v2, v3;
  List *flist, *nlist, *klist;
  double wgt[3][3];

  cmax = cost = 0.0;
  cmin = 1.0e20;
  count = 0;

  bool all_tet = true;
  for (i=0; i < node[n].element->max && all_tet; i++)
  {
    e = node[n].element->list[i];
    if (element_type(e) != 4)
      all_tet = false;
  }

  flist = new List();
  nlist = new List();
  klist = new List();

  // create list of nodes to test
  klist->Add_To_List(n);
  if (Exp_cost && !all_tet)
    for (i=0; i < node[n].nabor->max; i++)
      klist->Add_To_List(node[n].nabor->list[i]);

  for (h=0; h < klist->max; h++)
  {
    k = klist->list[h];
    for (i=0; i < node[n].element->max; i++)
    {
      e = node[n].element->list[i];

      if (adjacent_nodes(e,k,nlist,flist) <= 0)
        continue;

      jmax = (nlist->max == 3) ? 1 : nlist->max;

      for (j=0; j < jmax; j++)
      {
        f1 = flist->list[j];
        n1 = nlist->list[j];
        if (j < nlist->max-1)
        {
          f2 = flist->list[j+1];
          n2 = nlist->list[j+1];
        } else
        {
          f2 = flist->list[0];
          n2 = nlist->list[0];
        }
        if (j == 0)
        {
          f3 = flist->list[flist->max-1];
          n3 = nlist->list[nlist->max-1];
        } else
        {
          f3 = flist->list[j-1];
          n3 = nlist->list[j-1];
        }
        if (nlist->max > 3)
          f3 = 0;
        p0 = node[k].vert;
        p1 = node[n1].vert;
        p2 = node[n2].vert;
        p3 = node[n3].vert;

        v1 = Vector(p0,p1);
        v2 = Vector(p0,p2);
        v3 = Vector(p0,p3);
        J = scalar_triple_product(v1,v2,v3);
        if (J <= Jzero)
          c = MAX(1.0 - J,1.0+1.0e-12);
        else
        {
          double cn=1.0;
          if (all_tet)
            cn = tetrahedral_condition_number(p0,p1,p2,p3,1);
          else
          {
            double a1p, a1i, a2p, a2i, a3p, a3i;
            double d1p, d1i, d2p, d2i, d3p, d3i;
            d1p = d1i = d2p = d2i = d3p = d3i = 1.0;
            a1p = a1i = a2p = a2i = a3p = a3i = pi/3.0;
            a1i = pi - 2*pi/face[abs(f1)].node_list->max;
            a2i = pi - 2*pi/face[abs(f2)].node_list->max;
            if (f3 == 0)
              a3i = pi/3.0;
            else
              a3i = pi - 2*pi/face[abs(f3)].node_list->max;
            double c1 = cos(a1i);
            double c2 = cos(a2i);
            double c3 = cos(a3i);
            double s1 = sin(a1i);
            wgt[0][0] = 1.0;
            wgt[1][0] = 0.0;
            wgt[2][0] = 0.0;
            wgt[0][1] = c1;
            wgt[1][1] = s1;
            wgt[2][1] = 0.0;
            wgt[0][2] = c3;
            wgt[1][2] = (c2-c1*c3)/s1;
            wgt[2][2] = sqrt(MAX(0.0,1.0-wgt[0][2]*wgt[0][2]-wgt[1][2]*wgt[1][2]));
            cn = condition_number(p0,p1,p2,p3,wgt);
          }
          cn = MAX(1.0,cn);

          c = MAX(0.0,1.0 - 1.0/cn);
        }
        cost += c;
        cmax = MAX(cmax,c);
        count++;
      }
    }
  }
  delete flist;
  delete nlist;
  delete klist;

  cost/=MAX(1,count);
  double ratio = pow(cmax,(double)blend_exponent);
  cost=ratio*cmax + MAX(0.0,(1.0-ratio))*cost;
  return (cost);
}
