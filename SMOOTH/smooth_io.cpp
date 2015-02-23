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
#include "smooth.h"
#include "CGNS.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/
extern int rotate;
extern int displace;

void create_hybrid_maps(int nn, int nb, int ntet, int npyr, int npri, int nhex, int nply,
                        int **tet_n, int **pyr_n, int **pri_n, int **hex_n, int ***poly_n,
                        int *nt, int ***t_n, int *nq, int ***q_n, int *ngon, int ***ngon_n,
                        int **nmap, int **tri_map, int **quad_map, int **ngon_map,
                        int *tet_map, int *pyr_map, int *pri_map, int *hex_map, int *poly_map)
{
  // THIS ASSUMES THE NODE OWNERSHIP IS CORRECT!!!
  int b, e, i, j, n, o;
#ifdef PARALLEL
  int k, l, m, p;
  MPI_Status single_status;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt;
  int *recvcnt;
  MPI_Request *srequest;
  MPI_Request *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff;
  char **rbuff;
  List **pelems;
#endif

  // first identify owning processors for each element as
  // lowest processor from any of the nodes of the element
  for (b=0; b < nb; b++)
  {
    for (i=0; i < nt[b]; i++)
    {
      o = num_procs;
      for (j=0; j < 3; j++)
      {
        n = t_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        tri_map[b][i] = i;
      else
        tri_map[b][i] = -(o+1);
    }
    for (i=0; i < nq[b]; i++)
    {
      o = num_procs;
      for (j=0; j < 4; j++)
      {
        n = q_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        quad_map[b][i] = i;
      else
        quad_map[b][i] = -(o+1);
    }
    for (i=0; i < ngon[b]; i++)
    {
      o = num_procs;
      for (j=1; j <= ngon_n[b][i][0]; j++)
      {
        n = ngon_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        ngon_map[b][i] = i;
      else
        ngon_map[b][i] = -(o+1);
    }
  }

  for (e=0; e < ntet; e++)
  {
    o = num_procs;
    for (i=0; i < 4; i++)
    {
      n = tet_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      tet_map[e] = e;
    else
      tet_map[e] = -(o+1);
  }
  for (e=0; e < npyr; e++)
  {
    o = num_procs;
    for (i=0; i < 5; i++)
    {
      n = pyr_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      pyr_map[e] = e;
    else
      pyr_map[e] = -(o+1);
  }
  for (e=0; e < npri; e++)
  {
    o = num_procs;
    for (i=0; i < 6; i++)
    {
      n = pri_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      pri_map[e] = e;
    else
      pri_map[e] = -(o+1);
  }
  for (e=0; e < nhex; e++)
  {
    o = num_procs;
    for (i=0; i < 8; i++)
    {
      n = hex_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      hex_map[e] = e;
    else
      hex_map[e] = -(o+1);
  }
  for (e=0; e < nply; e++)
  {
    o = num_procs;
    for (i=1; i <= poly_n[e][0][0]; i++)
      for (j=1; j <= poly_n[e][i][0]; j++)
        o = MIN(o,nmap[poly_n[e][i][j]][1]);
    if (o == my_rank)
      poly_map[e] = e;
    else
      poly_map[e] = -(o+1);
  }
  
#ifdef PARALLEL
  // now take turns and set global element numbers
  if (num_procs > 1)
  {
    sendcnt = new int[num_procs];
    recvcnt = new int[num_procs];
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses  = new MPI_Status[num_procs];
    // allocate space for message buffers
    bdim = new int[num_procs];
    for (p=0; p < num_procs; p++)
      bdim[p] = 0;
    sbuff = new char*[num_procs];
    rbuff = new char*[num_procs];
    for (p=0; p < num_procs; p++)
    {
      sbuff[p] = 0;
      rbuff[p] = 0;
    }
    pelems = new List*[num_procs];
    for (p=0; p < num_procs; p++)
      pelems[p] = new List();
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int tet_total, pyr_total, pri_total, hex_total, poly_total;
  int *tri_total, *quad_total, *ngon_total;
  tri_total = new int[nb];
  quad_total = new int[nb];
  ngon_total = new int[nb];

  // Assign global element numbers
  if (my_rank == 0)
  {
    for (b=0; b < nb; b++)
    {
      tri_total[b] = quad_total[b] = ngon_total[b] = 0;
      for (i=0; i < nt[b]; i++)
        if (tri_map[b][i] >= 0)
          tri_map[b][i] = tri_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&tri_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&tri_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
      for (i=0; i < nq[b]; i++)
        if (quad_map[b][i] >= 0)
          quad_map[b][i] = quad_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&quad_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&quad_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
      for (i=0; i < ngon[b]; i++)
        if (ngon_map[b][i] >= 0)
          ngon_map[b][i] = ngon_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&ngon_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&ngon_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
    }

    tet_total = pyr_total = pri_total = hex_total = poly_total = 0;
    for (e=0; e < ntet; e++)
      if (tet_map[e] >= 0)
        tet_map[e] = tet_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&tet_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&tet_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of tet elements = %d",tet_total);
    //fflush(out_f);
    for (e=0; e < npyr; e++)
      if (pyr_map[e] >= 0)
        pyr_map[e] = pyr_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&pyr_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&pyr_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of pyr elements = %d",pyr_total);
    //fflush(out_f);
    for (e=0; e < npri; e++)
      if (pri_map[e] >= 0)
        pri_map[e] = pri_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&pri_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&pri_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of pri elements = %d",pri_total);
    //fflush(out_f);
    for (e=0; e < nhex; e++)
      if (hex_map[e] >= 0)
        hex_map[e] = hex_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&hex_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&hex_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of hex elements = %d",hex_total);
    //fflush(out_f);
    for (e=0; e < nply; e++)
      if (poly_map[e] >= 0)
        poly_map[e] = poly_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&poly_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&poly_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of poly elements = %d",poly_total);
    //fflush(out_f);
  } else
  {
#ifdef PARALLEL
    for (b=0; b < nb; b++)
    {
      MPI_Recv(&tri_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < nt[b]; i++)
        if (tri_map[b][i] >= 0)
          tri_map[b][i] = tri_total[b]++;
      MPI_Send(&tri_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
      MPI_Recv(&quad_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < nq[b]; i++)
        if (quad_map[b][i] >= 0)
          quad_map[b][i] = quad_total[b]++;
      MPI_Send(&quad_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
      MPI_Recv(&ngon_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < ngon[b]; i++)
        if (ngon_map[b][i] >= 0)
          ngon_map[b][i] = ngon_total[b]++;
      MPI_Send(&ngon_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    }
    MPI_Recv(&tet_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < ntet; e++)
      if (tet_map[e] >= 0)
        tet_map[e] = tet_total++;
    MPI_Send(&tet_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&pyr_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < npyr; e++)
      if (pyr_map[e] >= 0)
        pyr_map[e] = pyr_total++;
    MPI_Send(&pyr_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&pri_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < npri; e++)
      if (pri_map[e] >= 0)
        pri_map[e] = pri_total++;
    MPI_Send(&pri_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&hex_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < nhex; e++)
      if (hex_map[e] >= 0)
        hex_map[e] = hex_total++;
    MPI_Send(&hex_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&poly_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < nply; e++)
      if (poly_map[e] >= 0)
        poly_map[e] = poly_total++;
    MPI_Send(&poly_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
#endif
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);

  if (num_procs > 1)
  {
    // communicate and exchange the global number and indices not owned by me
    List **elist;
    int ne, ln;

    elist = new List*[nn];

    for (n=0; n < nn; n++)
      elist[n] = new List();

    int type, num_nodes;
    for (b=0; b < nb; b++)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      for (type=0; type < 3; type++)
      {
        for (n=0; n < nn; n++)
          elist[n]->Redimension(0);
        for (p=0; p < num_procs; p++)
          pelems[p]->Redimension(0);

        switch (type)
        {
          case 0:
            num_nodes = 3;
            for (e=0; e < nt[b]; e++)
            {
              for (i=0; i < num_nodes; i++)
              {
                n = t_n[b][e][i];
                elist[n]->Add_To_List(e);
              }
              if ((p=-tri_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
          case 1:
            num_nodes = 4;
            for (e=0; e < nq[b]; e++)
            {
              for (i=0; i < num_nodes; i++)
              {
                n = q_n[b][e][i];
                elist[n]->Add_To_List(e);
              }
              if ((p=-quad_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
          case 4:
            num_nodes = 0;
            for (e=0; e < ngon[b]; e++)
            {
              for (j=1; j <= ngon_n[b][e][0]; j++)
              {
                n = ngon_n[b][e][j];
                elist[n]->Add_To_List(e);
              }
              if ((p=-ngon_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
        }

        List nodelist;
        // send buffer size
        for (p=0; p < num_procs; p++)
        {
          sendcnt[p] = 0;
          if (p != my_rank)
          {
            sendcnt[p] = sizeof(int);
            for (i=0; i < pelems[p]->max; i++)
            {
              e = pelems[p]->list[i];
              switch (type)
              {
                case 0:
                  num_nodes = 3;
                  break;
                case 1:
                  num_nodes = 4;
                  break;
                case 2:
                  nodelist.Redimension(0);
                  for (j=1; j <= ngon_n[b][e][0]; j++)
                    nodelist.Check_List(ngon_n[b][e][j]);
                  num_nodes = nodelist.max;
                  break;
              }
              sendcnt[p] += (num_nodes*3 + 2)*sizeof(int);
            }
          }
        }
    
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          recvcnt[p] = 0;
          if (p != my_rank)
          {
            MPI_Isend(&sendcnt[p],1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_s++;
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // size message buffers for each processor
        for (p=0; p < num_procs; p++)
        {
          bdim[p] = 0;
          if (p == my_rank || (sendcnt[p] == 0 && recvcnt[p] == 0))
            continue;

          bdim[p] = MAX(sendcnt[p],recvcnt[p]);
          if (sbuff[p] > 0) delete[] sbuff[p];
          if (rbuff[p] > 0) delete[] rbuff[p];
          if (bdim[p] > 0) sbuff[p] = new char[bdim[p]];
          if (bdim[p] > 0) rbuff[p] = new char[bdim[p]];
        }

        // package the list of element per processor
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank || sendcnt[p] == 0)
            continue;

          sposition=0;
          MPI_Pack(&(pelems[p]->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          for (i=0; i < pelems[p]->max; i++)
          {
            e = pelems[p]->list[i];
            MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            switch(type)
            {
              case 0:
                num_nodes = 3;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < num_nodes; j++)
                {
                  n = t_n[b][e][j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
              case 1:
                num_nodes = 4;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < num_nodes; j++)
                {
                  n = q_n[b][e][j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
              case 2:
                nodelist.Redimension(0);
                for (j=1; j <= ngon_n[b][e][0]; j++)
                  nodelist.Check_List(ngon_n[b][e][j]);
                num_nodes = nodelist.max;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < nodelist.max; j++)
                {
                  n = nodelist.list[j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
            }
          }
        }

        // exchange messages
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
          if (sendcnt[p] > 0)
          {
            MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
          if (recvcnt[p] > 0)
          {
            MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // unpack element list and re-package the current map for requested elements
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;

          sposition=rposition=0;
          Pmap *lmap;
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Pack(&ne,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          for (i=0; i < ne; i++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&ln,1,MPI_INT,MPI_COMM_WORLD);
            lmap = new Pmap[ln];
            m = -1;
            for (n=0; n < ln; n++)
            {
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].global,1,MPI_INT,MPI_COMM_WORLD);
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].proc,1,MPI_INT,MPI_COMM_WORLD);
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].index,1,MPI_INT,MPI_COMM_WORLD);
              // find a node to use with the local node-element hash table
              if (m < 0 && lmap[n].proc == my_rank)
                m = lmap[n].index;
            }
            if (m < 0)
            {
              fprintf(out_f,"\nCreate_Hybrid_Maps: Boundary elements: no local node identified!");
              fflush(out_f);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            int le = -1; // local element identity
            for (j=0; j < elist[m]->max && le < 0; j++)
            {
              e = elist[m]->list[j];
              if (type == 2)
              {
                nodelist.Redimension(0);
                for (k=1; k <= ngon_n[b][e][0]; k++)
                  nodelist.Check_List(ngon_n[b][e][k]);
                if (nodelist.max != ln)
                  continue;
              }
              bool ematch = true;
              for (k=0; k < ln && ematch; k++)
              {
                switch(type)
                {
                  case 0: n = t_n[b][e][k]; break;
                  case 1: n = q_n[b][e][k]; break;
                  case 2: n = nodelist.list[k]; break;
                }
                bool nmatch = false;
                for (l=0; l < ln && !nmatch; l++)
                  if (lmap[l].proc == nmap[n][1] && lmap[l].index == nmap[n][2])
                    nmatch = true;
                ematch = nmatch;
              }
              if (ematch)
                le = e;
            }
            if (type == 2) nodelist.Redimension(0);

            if (le < 0)
            {
              fprintf(out_f,"\nCreate_Hybrid_Maps: Boundary Elements: no local element identified!");
              fflush(out_f);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }

            switch(type)
            {
              case 0: MPI_Pack(&tri_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
              case 1: MPI_Pack(&quad_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
              case 2: MPI_Pack(&ngon_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            }
  
            delete[] lmap;
          }
        }

        // exchange messages
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
          if (recvcnt[p] > 0)
          {
            MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
          if (sendcnt[p] > 0)
          {
            MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // unpack returned messages
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;

          sposition=rposition=0;
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
          for (m=0; m < ne; m++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
            switch(type)
            {
              case 0: tri_map[b][e] = i; break;
              case 1: quad_map[b][e] = i; break;
              case 2: ngon_map[b][e] = i; break;
            }
          }
        }
      }
    }

    for (type=0; type < 5; type++)
    {

      MPI_Barrier(MPI_COMM_WORLD);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: type = %d",type);
      //fflush(out_f);

      for (n=0; n < nn; n++)
        elist[n]->Redimension(0);
      for (p=0; p < num_procs; p++)
        pelems[p]->Redimension(0);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: elist and pelems reset");
      //fflush(out_f);

      // create node to element hash table
      // create list of elements to send to each processor
      switch (type)
      {
        case 0:
          num_nodes = 4;
          for (e=0; e < ntet; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = tet_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-tet_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 1:
          num_nodes = 5;
          for (e=0; e < npyr; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = pyr_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-pyr_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 2:
          num_nodes = 6;
          for (e=0; e < npri; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = pri_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-pri_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 3:
          num_nodes = 8;
          for (e=0; e < nhex; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = hex_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-hex_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 4:
          num_nodes = 0;
          for (e=0; e < nply; e++)
          {
            for (i=1; i <= poly_n[e][0][0]; i++)
              for (j=1; j <= poly_n[e][i][0]; j++)
              {
                n = poly_n[e][i][j];
                elist[n]->Add_To_List(e);
              }
            if ((p=-poly_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: pelems list made");
      //fflush(out_f);

      List nodelist;
      // send buffer size
      for (p=0; p < num_procs; p++)
      {
        sendcnt[p] = 0;
        if (p != my_rank)
        {
          sendcnt[p] = sizeof(int);
          for (i=0; i < pelems[p]->max; i++)
          {
            e = pelems[p]->list[i];
            switch (type)
            {
              case 0:
                num_nodes = 4;
                break;
              case 1:
                num_nodes = 5;
                break;
              case 2:
                num_nodes = 6;
                break;
              case 3:
                num_nodes = 8;
                break;
              case 4:
                nodelist.Redimension(0);
                for (j=1; j <= poly_n[e][0][0]; j++)
                  for (k=1; k <= poly_n[e][j][0]; k++)
                    nodelist.Check_List(poly_n[e][j][k]);
                num_nodes = nodelist.max;
                break;
            }
            sendcnt[p] += (num_nodes*3 + 2)*sizeof(int);
          }
        }
      }
    
      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: sendcnt made");
      //fflush(out_f);

      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        recvcnt[p] = 0;
        if (p != my_rank)
        {
          MPI_Isend(&sendcnt[p],1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_s++;
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: sendcnt exchanged");
      //fflush(out_f);

      // size message buffers for each processor
      for (p=0; p < num_procs; p++)
      {
        bdim[p] = 0;
        if (p == my_rank || (sendcnt[p] == 0 && recvcnt[p] == 0))
          continue;

        bdim[p] = MAX(sendcnt[p],recvcnt[p]);
        if (sbuff[p] > 0) delete[] sbuff[p];
        if (rbuff[p] > 0) delete[] rbuff[p];
        if (bdim[p] > 0) sbuff[p] = new char[bdim[p]];
        if (bdim[p] > 0) rbuff[p] = new char[bdim[p]];
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: buffers resized");
      //fflush(out_f);

      // package the list of element per processor
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank || sendcnt[p] == 0)
          continue;

        sposition=0;
        MPI_Pack(&(pelems[p]->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        for (i=0; i < pelems[p]->max; i++)
        {
          e = pelems[p]->list[i];
          MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          switch(type)
          {
            case 0:
              num_nodes = 4;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = tet_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 1:
              num_nodes = 5;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = pyr_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 2:
              num_nodes = 6;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = pri_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 3:
              num_nodes = 8;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = hex_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 4:
              nodelist.Redimension(0);
              for (j=1; j <= poly_n[e][0][0]; j++)
                for (k=1; k <= poly_n[e][j][0]; k++)
                  nodelist.Check_List(poly_n[e][j][k]);
              num_nodes = nodelist.max;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < nodelist.max; j++)
              {
                n = nodelist.list[j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
          }
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages packed");
      //fflush(out_f);

      // exchange messages
      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;
        if (sendcnt[p] > 0)
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }
        if (recvcnt[p] > 0)
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages exhanged");
      //fflush(out_f);

      // unpack element list and re-package the current map for requested elements
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        sposition=rposition=0;
        Pmap *lmap;
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Pack(&ne,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        for (i=0; i < ne; i++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ln,1,MPI_INT,MPI_COMM_WORLD);
          lmap = new Pmap[ln];
          m = -1;
          for (n=0; n < ln; n++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].global,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].proc,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].index,1,MPI_INT,MPI_COMM_WORLD);
            // find a node to use with the local node-element hash table
            if (m < 0 && lmap[n].proc == my_rank)
              m = lmap[n].index;
          }
          if (m < 0)
          {
            fprintf(out_f,"\nCreate_Hybrid_Maps: no local node identified!");
            fflush(out_f);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }

          int le = -1; // local element identity
          for (j=0; j < elist[m]->max && le < 0; j++)
          {
            e = elist[m]->list[j];
            if (type == 4)
            {
              nodelist.Redimension(0);
              for (j=1; j <= poly_n[e][0][0]; j++)
                for (k=1; k <= poly_n[e][j][0]; k++)
                  nodelist.Check_List(poly_n[e][j][k]);
              if (nodelist.max != ln)
                continue;
            }
            bool ematch = true;
            for (k=0; k < ln && ematch; k++)
            {
              switch(type)
              {
                case 0: n = tet_n[e][k]; break;
                case 1: n = pyr_n[e][k]; break;
                case 2: n = pri_n[e][k]; break;
                case 3: n = hex_n[e][k]; break;
                case 4: n = nodelist.list[k]; break;
              }
              bool nmatch = false;
              for (l=0; l < ln && !nmatch; l++)
                if (lmap[l].proc == nmap[n][1] && lmap[l].index == nmap[n][2])
                  nmatch = true;
              ematch = nmatch;
            }
            if (ematch)
              le = e;
          }
          if (type == 4) nodelist.Redimension(0);

          if (le < 0)
          {
            fprintf(out_f,"\nCreate_Hybrid_Maps: no local element identified!");
            fflush(out_f);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }

          switch(type)
          {
            case 0: MPI_Pack(&tet_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 1: MPI_Pack(&pyr_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 2: MPI_Pack(&pri_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 3: MPI_Pack(&hex_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 4: MPI_Pack(&poly_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
          }

          delete[] lmap;
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages repacked");
      //fflush(out_f);

      // exchange messages
      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;
        if (recvcnt[p] > 0)
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }
        if (sendcnt[p] > 0)
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages re-exchanged");
      //fflush(out_f);

      // unpack returned messages
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        sposition=rposition=0;
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
        for (m=0; m < ne; m++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
          switch(type)
          {
            case 0: tet_map[e] = i; break;
            case 1: pyr_map[e] = i; break;
            case 2: pri_map[e] = i; break;
            case 3: hex_map[e] = i; break;
            case 4: poly_map[e] = i; break;
          }
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages unpacked and saved");
      //fflush(out_f);

    }

    // release communication memory
    for (p=0; p < num_procs; p++)
    {
      if (sbuff[p] > 0) delete[] sbuff[p];
      if (rbuff[p] > 0) delete[] rbuff[p];
      delete pelems[p];
    }
    delete[] pelems;
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] srequest;
    delete[] rrequest;
    delete[] statuses;
    // allocate space for message buffers
    delete[] bdim;
    delete[] sbuff;
    delete[] rbuff;
    for (n=0; n < nn; n++)
      delete elist[n];
    delete[] elist;
  }
#endif

  //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: done");
  //fflush(out_f);

  return;
}

void POLYMESH::exchange_node_int(int *array)
{
#ifdef PARALLEL
  int temp_node, temp_node2, n, p;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  
  sendcnt = new int [num_procs];
  recvcnt = new int [num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sendcnt[n] = 0;
    recvcnt[n] = 0;
  }
  srequest = new MPI_Request[num_procs];
  rrequest = new MPI_Request[num_procs];
  statuses = new MPI_Status[num_procs];

  //cycle thru nodes
  for (n = 0; n < nn; n++)
    if (node[n].me.proc != my_rank)
      sendcnt[node[n].me.proc]++; 
    
  //now, send and recv count for all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (p == my_rank)
      continue;

    MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
    nreq_s++;
 
    MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
    nreq_r++;
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //indices
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (node[n].me.proc == p)
      {
        MPI_Pack(&(node[n].me.index),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
 
  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  
  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      //pack node and nrm
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(array[temp_node]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(array[temp_node]),1,MPI_INT,MPI_COMM_WORLD);
    }
  }

  //finally, free MPI mem
  delete[] sendcnt;
  delete[] recvcnt;
  delete[] srequest;
  delete[] rrequest;
  delete[] statuses;
  for (n = 0; n < num_procs; n++)
  {
    if (sbuff[n] != 0) delete [] sbuff[n];
    if (rbuff[n] != 0) delete [] rbuff[n];
  }
  delete[] sbuff;
  delete[] rbuff;
  delete[] bdim;
#endif

  return;
}

void POLYMESH::exchange_node_double(double *array)
{
#ifdef PARALLEL
  int temp_node, temp_node2, n, p;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  
  sendcnt = new int [num_procs];
  recvcnt = new int [num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sendcnt[n] = 0;
    recvcnt[n] = 0;
  }
  srequest = new MPI_Request[num_procs];
  rrequest = new MPI_Request[num_procs];
  statuses = new MPI_Status[num_procs];

  //cycle thru nodes
  for (n = 0; n < nn; n++)
    if (node[n].me.proc != my_rank)
      sendcnt[node[n].me.proc]++; 
    
  //now, send and recv count for all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (p == my_rank)
      continue;

    MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
    nreq_s++;
 
    MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
    nreq_r++;
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*(sizeof(int)+sizeof(double)); //indices
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (node[n].me.proc == p)
      {
        MPI_Pack(&(node[n].me.index),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
 
  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  
  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      //pack node and nrm
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(array[temp_node]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(array[temp_node]),1,MPI_DOUBLE,MPI_COMM_WORLD);
    }
  }

  //finally, free MPI mem
  delete[] sendcnt;
  delete[] recvcnt;
  delete[] srequest;
  delete[] rrequest;
  delete[] statuses;
  for (n = 0; n < num_procs; n++)
  {
    if (sbuff[n] != 0) delete [] sbuff[n];
    if (rbuff[n] != 0) delete [] rbuff[n];
  }
  delete[] sbuff;
  delete[] rbuff;
  delete[] bdim;
#endif

  return;
}

void POLYMESH::exchange_node_vector(Vector *array)
{
#ifdef PARALLEL
  int temp_node, temp_node2, n, p;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff, **rbuff;
  
  sendcnt = new int [num_procs];
  recvcnt = new int [num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sendcnt[n] = 0;
    recvcnt[n] = 0;
  }
  srequest = new MPI_Request[num_procs];
  rrequest = new MPI_Request[num_procs];
  statuses = new MPI_Status[num_procs];

  //cycle thru nodes
  for (n = 0; n < nn; n++)
    if (node[n].me.proc != my_rank)
      sendcnt[node[n].me.proc]++; 
    
  //now, send and recv count for all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (p == my_rank)
      continue;

    MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
    nreq_s++;
 
    MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
    nreq_r++;
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //allocate buffers to send/recv node indices
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*(sizeof(int)+3*sizeof(double)); //indices
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }

  //now, pack node indices as per procs that will be sending them
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nn; n++)
    {
      if (node[n].me.proc == p)
      {
        MPI_Pack(&(node[n].me.index),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
 
  temp_node = 0; //to hold recv proc local index
  temp_node2 = 0; //to hold send proc local index
  
  //now, unpack local indices and repack nodes with other index
  for (p = 0; p < num_procs; p++)
  {
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD); 
      //unpack index to send back
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
      //pack node and nrm
      MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(array[temp_node][0]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(array[temp_node][1]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      MPI_Pack(&(array[temp_node][2]),1,MPI_DOUBLE,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
    }
  }

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack index and node
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(array[temp_node][0]),1,MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(array[temp_node][1]),1,MPI_DOUBLE,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&(array[temp_node][2]),1,MPI_DOUBLE,MPI_COMM_WORLD);
    }
  }

  //finally, free MPI mem
  delete[] sendcnt;
  delete[] recvcnt;
  delete[] srequest;
  delete[] rrequest;
  delete[] statuses;
  for (n = 0; n < num_procs; n++)
  {
    if (sbuff[n] != 0) delete [] sbuff[n];
    if (rbuff[n] != 0) delete [] rbuff[n];
  }
  delete[] sbuff;
  delete[] rbuff;
  delete[] bdim;
#endif

  return;
}

int POLYMESH::smooth_io(int mode, geometry *geom, char sname[])
{
  int flag;
  
  flag = 0;

  switch (mode)
  {
    case -1:
      //read in comp and physical meshes
      if (read_mesh(sname) != 0)
      {
        fprintf(out_f,"\nCouldn't open file %s",sname);
        fflush(out_f);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }

      if (nb != geom->ngb)
      {
        fprintf(out_f,"\nNumber of boundaries in mesh does not match geometry file!");
        fflush(out_f);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
      break;

    case 1:

      // write new physical mesh file

      //mesh_stats(tot_neg,tot_neg_skw,tot_pos_skw,tot_pos);
  
      if (write_mesh(sname) != 0)
      {
        fprintf(out_f,"\nCouldn't open file %s",sname);
        fflush(out_f);
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD,0);
#endif
        exit(0);
      }
      break;
  }

  return flag;
}

int POLYMESH::face_number(List *nlist, List **flist)
{
  int i, j, m, n0, n1, f, fnum;
  bool test;

  if (flist == 0 || nlist == 0)
  {
    fprintf(stderr,"\nFACE_NUMBER: node list and/or face list not defined!");
    exit(0);
  }

  fnum = 0;
  test = true;
  if (nlist->max == 0) test = false;
  for (i=0; i < nlist->max && test; i++)
  {
    m = nlist->list[i];
    if (flist[m]->max == 0) test = false;
  }
  if (test)
  {
    n0 = nlist->list[0];
    for (i=0; i < flist[n0]->max && fnum==0; i++)
    {
      f = flist[n0]->list[i];
      if (f <= 0)
      {
        fprintf(stderr,"\nFACE_NUMBER: Zero or negative face number stored in flist.");
        exit(0);
      }
      test = true;
      for (j=1; j < nlist->max && test; j++)
      {
        n1 = nlist->list[j];
        if (!flist[n1]->Is_In_List(f))
          test = false;
      }
      if (test)
        fnum = f;
    }
  }

  if (fnum > 0 && fnum <= nfaces && nlist->max == face[fnum].node_list->max)
  {
    int wound = 0;
    for (i=0; i < face[fnum].node_list->max && wound == 0; i++)
    {
      n0 = face[fnum].node_list->list[i];
      n1 = face[fnum].node_list->list[(i+1)%face[fnum].node_list->max];
      if (n0 == nlist->list[0] && n1 == nlist->list[1])
        wound = 1;
      else if (n1 == nlist->list[0] && n0 == nlist->list[1])
        wound = -1;
    }
    if (wound == 0)
    {
      fprintf(stderr,"\nFACE_NUMBER: Unable to determine winding!");
      exit(0);
    }
    fnum *=wound;
  } else if (fnum == 0)
  {
    if (nfaces+1 >= face_dim)
    {
      int chunk = 100000;
      face_dim += chunk;
      if (num_procs == 1)
        fprintf(out_f,"\nFACE_NUMBER: Face dimension about to be exceeded, growing by %d to %d",chunk,face_dim);
      face = (POLYGONAL_ELEMENT*)realloc((void*)face,face_dim*sizeof(POLYGONAL_ELEMENT));
      for (i=nfaces+1; i < face_dim; i++)
        face[i].Initialize();
    }
    face[fnum=++nfaces].Add_Nodes(nlist);
    for (m=0; m < nlist->max; m++)
      flist[nlist->list[m]]->Add_To_List(fnum);
  }
  return(fnum);
}

int POLYMESH::read_mesh(char sname[]) //open and read mesh CGNS file
{
  int b, n, i, j, k;
  int flag = 0;
  int parallel = 1;

  // THIS ROUTINE CURRENTLY READS FOUR BASIC ELEMENTS AND CONVERTS TO POLYHEDRA
  // NEED TO ADD POLYHEDRA ELEMENTS TO P_CGNS ROUTINE!!!!!

  int ntet, npyr, npri, nhex, nply;
  Point *p;
  char **bname;
  int *nt, ***t_n, *nq, ***q_n, *ng, ***ngon_n;
  int **tet_n, **pyr_n, **pri_n, **hex_n, ***poly_n;
  int **nmap, **tri_map, **quad_map, **ngon_map;
  int *tet_map, *pyr_map, *pri_map, *hex_map, *poly_map;
  int g_nn, g_nfaces, g_nelems, g_ntet, g_npyr, g_npri, g_nhex, g_nply, *g_nt, *g_nq, *g_ng; // global counts (excluding ghosts)

  nn = nb = ntet = npyr = npri = nhex = nply = -1;
  nt = 0;
  t_n = 0;
  nq = 0;
  q_n = 0;
  ng = 0;
  ngon_n = 0;
  tet_n = 0;
  pyr_n = 0;
  pri_n = 0;
  hex_n = 0;
  poly_n = 0;
  nmap = 0;
  tri_map = 0;
  quad_map = 0;
  ngon_map = 0;
  tet_map = 0;
  pyr_map = 0;
  pri_map = 0;
  hex_map = 0;
  poly_map = 0;
  p = 0;

  //actually reads in files and allocates memory for the object
  if (num_procs > 1)
  {
    flag = P_CGNS_read(parallel,sname,nn,&p,nb,&bname,&nt,&t_n,&nq,&q_n,&ng,&ngon_n,ntet,&tet_n,
                       npyr,&pyr_n,npri,&pri_n,nhex,&hex_n,nply,&poly_n,&nmap,&tri_map,&quad_map,&ngon_map,
                       &tet_map,&pyr_map,&pri_map,&hex_map,&poly_map);
  } else
  {
    //no need for maps since no need to reassemble
    flag = CGNS_read(sname,nn,&p,nb,&bname,&nt,&t_n,&nq,&q_n,&ng,&ngon_n,ntet,&tet_n,npyr,&pyr_n,
                     npri,&pri_n,nhex,&hex_n,nply,&poly_n);
  }

  //fprintf(out_f,"\nProcess %d: Local number of nodes      = %d",my_rank,nn);
  //fflush(out_f);
  //if (ntet > 0) fprintf(out_f,"\nProcess %d: Local number of tetrahedra = %d",my_rank,ntet);
  //fflush(out_f);
  //if (npyr > 0) fprintf(out_f,"\nProcess %d: Local number of pyramid    = %d",my_rank,npyr);
  //fflush(out_f);
  //if (npri > 0) fprintf(out_f,"\nProcess %d: Local number of prism      = %d",my_rank,npri);
  //fflush(out_f);
  //if (nhex > 0) fprintf(out_f,"\nProcess %d: Local number of hexahedra  = %d",my_rank,nhex);
  //fflush(out_f);
  //if (nply > 0) fprintf(out_f,"\nProcess %d: Local number of polyhedra  = %d",my_rank,nply);
  //fflush(out_f);
  //fprintf(out_f,"\nProcess %d: Local number of boundaries = %d",my_rank,nb);
  //fflush(out_f);

  nelems = ntet+npyr+npri+nhex+nply;
  nfaces = ntet*4+npyr*5+npri*5+nhex*6;
  for (n=0; n < nply; n++)
    nfaces += poly_n[n][0][0];
  
  for (n=0; n < nb; n++)
  {
    nfaces += (nt[n]+nq[n]+ng[n]);
    //if (nt[n] > 0)
    //  fprintf(out_f,"\nProcess %d: Boundary %d, Local # of triangles      = %d",my_rank,n+1,nt[n]);
    //fflush(out_f);
    //if (nq[n] > 0)
    //  fprintf(out_f,"\nProcess %d: Boundary %d, Local # of quadrilaterals = %d",my_rank,n+1,nq[n]);
    //fflush(out_f);
    //if (ng[n] > 0)
    //  fprintf(out_f,"\nProcess %d: Boundary %d, Local # of polygons       = %d",my_rank,n+1,ng[n]);
    //fflush(out_f);
  }
  nfaces /= 2;
  
  //fprintf(out_f,"\n");
  //fflush(out_f);
  //fprintf(out_f,"\nProcess %d: Local number of nodes      = %d",my_rank,nn);
  //fflush(out_f);
  //fprintf(out_f,"\nProcess %d: Local number of elements   = %d",my_rank,nelems);
  //fflush(out_f);

  node_dim = nn;
  element_dim = nelems;
  face_dim = nfaces+1; //FACES ARE NUMBERED FROM 1 not 0 !!!!!
  boundary_dim = nb;
  node = (SNODE*)malloc(node_dim*sizeof(SNODE));
  element = (POLYHEDRAL_ELEMENT*)malloc(element_dim*sizeof(POLYHEDRAL_ELEMENT));
  face = (POLYGONAL_ELEMENT*)malloc(face_dim*sizeof(POLYGONAL_ELEMENT));
  boundary = (BOUNDARY_MESH*)malloc(boundary_dim*sizeof(BOUNDARY_MESH));
  for (i=0; i < element_dim; i++)
    element[i].Initialize();
  for (i=0; i < face_dim; i++)
    face[i].Initialize();
  for (i=0; i < boundary_dim; i++)
    boundary[i].Initialize();
  
#ifdef PARALLEL
  // determine global numbers
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  g_nn = 0;
  int nmax = 0;
  for (i = 0; i < nn; i++)
  {
    node[i].vert = p[i];
    if (num_procs > 1)
    {
      node[i].me = Pmap(nmap[i][0],nmap[i][1],nmap[i][2]);
      free(nmap[i]);
    } else
    {
      node[i].me = Pmap(i,0,i);
    }
    nmax = MAX(nmax,node[i].me.global+1);
  }
  free(p);
  if (num_procs > 1)
    free(nmap);

  g_nn=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_nn,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

g_nfaces = nfaces;
#ifdef PARALLEL
  MPI_Allreduce(&nfaces,&g_nfaces,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

  nmax = 0;
  if (num_procs == 1)
    nmax = ntet;
  else
    for (i=0; i < ntet; i++)
      nmax = MAX(nmax,tet_map[i]+1);

  g_ntet=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_ntet,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  nmax = 0;
  if (num_procs == 1)
    nmax = npyr;
  else
    for (i=0; i < npyr; i++)
      nmax = MAX(nmax,pyr_map[i]+1);

  g_npyr=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_npyr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  nmax = 0;
  if (num_procs == 1)
    nmax = npri;
  else
    for (i=0; i < npri; i++)
      nmax = MAX(nmax,pri_map[i]+1);

  g_npri=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_npri,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  nmax = 0;
  if (num_procs == 1)
    nmax = nhex;
  else
    for (i=0; i < nhex; i++)
      nmax = MAX(nmax,hex_map[i]+1);

  g_nhex=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_nhex,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  nmax = 0;
  if (num_procs == 1)
    nmax = nply;
  else
    for (i=0; i < nply; i++)
      nmax = MAX(nmax,poly_map[i]+1);

  g_nply=nmax;
#ifdef PARALLEL
  MPI_Allreduce(&nmax,&g_nply,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  g_nt = (int*)malloc(nb*sizeof(int));
  g_nq = (int*)malloc(nb*sizeof(int));
  g_ng = (int*)malloc(nb*sizeof(int));

  for (b=0; b < nb; b++)
  {
    g_nt[b]=nt[b];
    g_nq[b]=nq[b];
    g_ng[b]=ng[b];
#ifdef PARALLEL
    MPI_Allreduce(&nt[b],&g_nt[b],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nq[b],&g_nq[b],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&ng[b],&g_ng[b],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  }
  
  if (my_rank == 0)
  {
    fprintf(out_f,"\n\nGlobal number of nodes      = %d",g_nn);
    fprintf(out_f,"\n\nGlobal number of faces      = %d",g_nfaces);
    if (g_ntet > 0) fprintf(out_f,"\nGlobal number of tetrahedra = %d",g_ntet);
    if (g_npyr > 0) fprintf(out_f,"\nGlobal number of pyramid    = %d",g_npyr);
    if (g_npri > 0) fprintf(out_f,"\nGlobal number of prism      = %d",g_npri);
    if (g_nhex > 0) fprintf(out_f,"\nGlobal number of hexahedra  = %d",g_nhex);
    if (g_nply > 0) fprintf(out_f,"\nGlobal number of polyhedra  = %d",g_nply);
    fprintf(out_f,"\nGlobal number of boundaries = %d",nb);

    for (n=0; n < nb; n++)
    {
      fprintf(out_f,"\nBoundary %d:",n+1);
      if (g_nt[n] > 0)
        fprintf(out_f,"\n Global # of triangles      = %d",g_nt[n]);
      if (g_nq[n] > 0)
        fprintf(out_f,"\n Global # of quadrilaterals = %d",g_nq[n]);
      if (g_ng[n] > 0)
        fprintf(out_f,"\n Global # of polygons       = %d",g_ng[n]);
    }
    fprintf(out_f,"\n");
    fflush(out_f);
  }

  List *lface, **flist;
  flist = new List*[nn];
  for (n=0; n < nn; n++)
    flist[n] = new List();

  lface = new List();
  nelems = 0;
  int fnum;
  nfaces = 0;

  g_nelems = 0;
  // add tetrahdedra to element array
  lface->Redimension(3);
  lface->max = 3;
  for (i=0; i < ntet; i++)
  {
    element[nelems].Initialize();
    element[nelems].node_list->Redimension(4);
    element[nelems].node_list->Add_To_List(tet_n[i][0]);
    element[nelems].node_list->Add_To_List(tet_n[i][1]);
    element[nelems].node_list->Add_To_List(tet_n[i][2]);
    element[nelems].node_list->Add_To_List(tet_n[i][3]);

    lface->list[0]=tet_n[i][0];
    lface->list[1]=tet_n[i][2];
    lface->list[2]=tet_n[i][1];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=tet_n[i][0];
    lface->list[1]=tet_n[i][1];
    lface->list[2]=tet_n[i][3];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=tet_n[i][1];
    lface->list[1]=tet_n[i][2];
    lface->list[2]=tet_n[i][3];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=tet_n[i][2];
    lface->list[1]=tet_n[i][0];
    lface->list[2]=tet_n[i][3];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    if (num_procs > 1)
      element[nelems].me = Pmap(tet_map[i] + g_nelems,my_rank,nelems);  // NOT CORRECT!!
    else
      element[nelems].me = Pmap(nelems,0,nelems);
    nelems++;
  }
  g_nelems += g_ntet;
  // free tetrahedra memory
  if (ntet > 0)
  {
    for (i=0; i < ntet; i++)
      free(tet_n[i]);
    free(tet_n);
    if (num_procs > 1)
      free(tet_map);
    tet_n=0;
    ntet=0;
  }
    
  // add pyramid to element array
  lface->Redimension(4);
  for (i=0; i < npyr; i++)
  {
    element[nelems].Initialize();
    element[nelems].node_list->Redimension(5);
    element[nelems].node_list->Add_To_List(pyr_n[i][0]);
    element[nelems].node_list->Add_To_List(pyr_n[i][1]);
    element[nelems].node_list->Add_To_List(pyr_n[i][2]);
    element[nelems].node_list->Add_To_List(pyr_n[i][3]);
    element[nelems].node_list->Add_To_List(pyr_n[i][4]);

    lface->max = 4;
    lface->list[0]=pyr_n[i][0];
    lface->list[1]=pyr_n[i][3];
    lface->list[2]=pyr_n[i][2];
    lface->list[3]=pyr_n[i][1];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->max = 3;
    lface->list[0]=pyr_n[i][0];
    lface->list[1]=pyr_n[i][1];
    lface->list[2]=pyr_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pyr_n[i][1];
    lface->list[1]=pyr_n[i][2];
    lface->list[2]=pyr_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pyr_n[i][2];
    lface->list[1]=pyr_n[i][3];
    lface->list[2]=pyr_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pyr_n[i][3];
    lface->list[1]=pyr_n[i][0];
    lface->list[2]=pyr_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    if (num_procs > 1)
      element[nelems].me = Pmap(pyr_map[i] + g_nelems,my_rank,nelems);  // NOT CORRECT!!
    else
      element[nelems].me = Pmap(nelems,0,nelems);
    nelems++;
  }
  g_nelems += g_npyr;
  // free pyramid memory
  if (npyr > 0)
  {
    for (i=0; i < npyr; i++)
      free(pyr_n[i]);
    free(pyr_n);
    if (num_procs > 1)
      free(pyr_map);
    pyr_n=0;
    npyr=0;
  }
    
  // add prisms to element array
  lface->Redimension(4);
  for (i=0; i < npri; i++)
  {
    element[nelems].Initialize();
    element[nelems].node_list->Redimension(6);
    element[nelems].node_list->Add_To_List(pri_n[i][0]);
    element[nelems].node_list->Add_To_List(pri_n[i][1]);
    element[nelems].node_list->Add_To_List(pri_n[i][2]);
    element[nelems].node_list->Add_To_List(pri_n[i][3]);
    element[nelems].node_list->Add_To_List(pri_n[i][4]);
    element[nelems].node_list->Add_To_List(pri_n[i][5]);

    lface->max = 4;
    lface->list[0]=pri_n[i][0];
    lface->list[1]=pri_n[i][1];
    lface->list[2]=pri_n[i][4];
    lface->list[3]=pri_n[i][3];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pri_n[i][1];
    lface->list[1]=pri_n[i][2];
    lface->list[2]=pri_n[i][5];
    lface->list[3]=pri_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pri_n[i][2];
    lface->list[1]=pri_n[i][0];
    lface->list[2]=pri_n[i][3];
    lface->list[3]=pri_n[i][5];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->max = 3;
    lface->list[0]=pri_n[i][0];
    lface->list[1]=pri_n[i][2];
    lface->list[2]=pri_n[i][1];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=pri_n[i][3];
    lface->list[1]=pri_n[i][4];
    lface->list[2]=pri_n[i][5];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    if (num_procs > 1)
      element[nelems].me = Pmap(pri_map[i] + g_nelems,my_rank,nelems);  // NOT CORRECT!!
    else
      element[nelems].me = Pmap(nelems,0,nelems);
    nelems++;
  }
  g_nelems += g_npri;
  // free prism memory
  if (npri > 0)
  {
    for (i=0; i < npri; i++)
      free(pri_n[i]);
    free(pri_n);
    if (num_procs > 1)
      free(pri_map);
    pri_n=0;
    npri=0;
  }

  // add hexes to element array
  lface->Redimension(4);
  lface->max = 4;
  for (i=0; i < nhex; i++)
  {
    element[nelems].Initialize();
    element[nelems].node_list->Redimension(8);
    element[nelems].node_list->Add_To_List(hex_n[i][0]);
    element[nelems].node_list->Add_To_List(hex_n[i][1]);
    element[nelems].node_list->Add_To_List(hex_n[i][2]);
    element[nelems].node_list->Add_To_List(hex_n[i][3]);
    element[nelems].node_list->Add_To_List(hex_n[i][4]);
    element[nelems].node_list->Add_To_List(hex_n[i][5]);
    element[nelems].node_list->Add_To_List(hex_n[i][6]);
    element[nelems].node_list->Add_To_List(hex_n[i][7]);

    lface->list[0]=hex_n[i][0];
    lface->list[1]=hex_n[i][3];
    lface->list[2]=hex_n[i][2];
    lface->list[3]=hex_n[i][1];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=hex_n[i][0];
    lface->list[1]=hex_n[i][1];
    lface->list[2]=hex_n[i][5];
    lface->list[3]=hex_n[i][4];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=hex_n[i][1];
    lface->list[1]=hex_n[i][2];
    lface->list[2]=hex_n[i][6];
    lface->list[3]=hex_n[i][5];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=hex_n[i][2];
    lface->list[1]=hex_n[i][3];
    lface->list[2]=hex_n[i][7];
    lface->list[3]=hex_n[i][6];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=hex_n[i][0];
    lface->list[1]=hex_n[i][4];
    lface->list[2]=hex_n[i][7];
    lface->list[3]=hex_n[i][3];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    lface->list[0]=hex_n[i][4];
    lface->list[1]=hex_n[i][5];
    lface->list[2]=hex_n[i][6];
    lface->list[3]=hex_n[i][7];
    fnum = face_number(lface,flist);
    element[nelems].Add_Face(fnum);

    if (num_procs > 1)
      element[nelems].me = Pmap(hex_map[i] + g_nelems,my_rank,nelems);  // NOT CORRECT!!
    else
      element[nelems].me = Pmap(nelems,0,nelems);
    nelems++;
  }
  g_nelems += g_nhex;
  // free hex memory
  if (nhex > 0)
  {
    for (i=0; i < nhex; i++)
      free(hex_n[i]);
    free(hex_n);
    if (num_procs > 1)
      free(hex_map);
    hex_n=0;
    nhex=0;
  }

  // add polyhedra to element array
  lface->Redimension(0);
  for (i=0; i < nply; i++)
  {
    element[nelems].Initialize();
    for (j=1; j <= poly_n[i][0][0]; j++)
    {
      lface->Redimension(poly_n[i][j][0]);
      lface->max = 0;
      for (k=1; k <= poly_n[i][j][0]; k++)
      {
        lface->Add_To_List(poly_n[i][j][k]);
        element[nelems].node_list->Check_List(poly_n[i][j][k]);
      }
      fnum = face_number(lface,flist);
      element[nelems].Add_Face(fnum);
    }
    if (num_procs > 1)
      element[nelems].me = Pmap(poly_map[i] + g_nelems,my_rank,nelems);  // NOT CORRECT!!
    else
      element[nelems].me = Pmap(nelems,0,nelems);
    element[nelems].node_list->Redimension(element[nelems].node_list->max);
    nelems++;
  }
  g_nelems += g_nply;
  if (nply > 0)
  {
    for (i=0; i < nply; i++)
    {
      for (j=poly_n[i][0][0]; j >= 0; j--)
        free(poly_n[i][j]);
      free(poly_n[i]);
    }
    poly_n=0;
    nply=0;
  }

  // process boundaries
  for (b=0; b < nb; b++)
  {
    boundary[b].Initialize();
    boundary[b].polygon_list = new List();
    boundary[b].name = (char*)malloc((strlen(bname[b])+1)*sizeof(char));
    strcpy(boundary[b].name,bname[b]);
    free(bname[b]);
    if (nt[b] > 0)
    {
      lface->Redimension(3);
      lface->max = 3;
      for (i=0; i < nt[b]; i++)
      {
        lface->list[0]=t_n[b][i][0];
        lface->list[1]=t_n[b][i][1];
        lface->list[2]=t_n[b][i][2];
        fnum = face_number(lface,flist);
        boundary[b].polygon_list->Add_To_List(fnum);
        free(t_n[b][i]);
      }
      free(t_n[b]);
      if (num_procs > 1)
        free(tri_map[b]);
      nt[b]=0;
    }
    if (nq[b] > 0)
    {
      lface->Redimension(4);
      lface->max = 4;
      for (i=0; i < nq[b]; i++)
      {
        lface->list[0]=q_n[b][i][0];
        lface->list[1]=q_n[b][i][1];
        lface->list[2]=q_n[b][i][2];
        lface->list[3]=q_n[b][i][3];
        fnum = face_number(lface,flist);
        boundary[b].polygon_list->Add_To_List(fnum);
        free(q_n[b][i]);
      }
      free(q_n[b]);
      if (num_procs > 1)
        free(quad_map[b]);
      nq[b]=0;
    }
    if (ng[b] > 0)
    {
      for (i=0; i < ng[b]; i++)
      {
        lface->Redimension(0);
        lface->max = 0;
        for (j=1; j <= ngon_n[b][i][0]; j++)
          lface->Add_To_List(ngon_n[b][i][j]);
        fnum = face_number(lface,flist);
        boundary[b].polygon_list->Add_To_List(fnum);
        free(ngon_n[b][i]);
      }
      free(ngon_n[b]);
      if (num_procs > 1)
        free(ngon_map[b]);
      ng[b]=0;
    }
  }
  free(g_nt);
  free(g_nq);
  free(g_ng);
  free(bname);
  free(nt);
  free(nq);
  free(ng);
  free(t_n);
  free(q_n);
  free(ngon_n);
  if (num_procs > 1)
  {
    free(tri_map);
    free(quad_map);
    free(ngon_map);
  }

  delete lface;
  for (n=0; n < nn; n++)
    delete flist[n];
  delete[] flist;

  create_element_maps();

  //return flag for success
  return (flag);
}

int POLYMESH::write_mesh(char sname[])
{
  int b, e, f, n, i, j, k;
  int flag = 0; 
  int parallel = 1;
  
  // THIS ROUTINE CURRENTLY WRITES FOUR BASIC ELEMENTS AND CONVERTS FROM POLYHEDRA

  int ntet, npyr, npri, nhex, nply;
  Point *p;
  char **bname;
  int *nt, ***t_n, *nq, ***q_n, *ng, ***ngon_n;
  int **tet_n, **pyr_n, **pri_n, **hex_n, ***poly_n;
  int **nmap, **tri_map, **quad_map, **ngon_map;
  int *tet_map, *pyr_map, *pri_map, *hex_map, *poly_map;

  ntet = npyr = npri = nhex = nply = 0;

  p = (Point*)malloc(nn*sizeof(Point));
  for (n=0; n < nn; n++)
    p[n] = node[n].vert;
  if (num_procs > 1)
  {
    nmap = (int**)malloc(nn*sizeof(int*));
    for (n=0; n < nn; n++)
    {
      nmap[n] = (int*)malloc(3*sizeof(int));
      nmap[n][0] = node[n].me.global;
      nmap[n][1] = node[n].me.proc;
      nmap[n][2] = node[n].me.index;
    }
  } else
    nmap = 0;

  nt = (int*)malloc(nb*sizeof(int));
  nq = (int*)malloc(nb*sizeof(int));
  ng = (int*)malloc(nb*sizeof(int));

  // count number of element of each type
  for (e=0; e < nelems; e++)
  {
    switch(element_type(e))
    {
      case 1: break;
      case 4: ntet++; break;
      case 5: npyr++; break;
      case 6: npri++; break;
      case 8: nhex++; break;
      default: nply++; break;
    }
  }
  bname = (char**)malloc(nb*sizeof(char*));
  for (b=0; b < nb; b++)
  {
    nt[b] = nq[b] = ng[b] = 0;
    bname[b] = (char*)malloc(33*sizeof(char));
    strcpy(bname[b],boundary[b].name);
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = abs(boundary[b].polygon_list->list[i]);
      if (face[f].node_list->max == 3)
        nt[b]++;
      if (face[f].node_list->max == 4)
        nq[b]++;
      if (face[f].node_list->max > 4)
        ng[b]++;
    }
  }
  
  tri_map=0;
  quad_map=0;
  ngon_map=0;
  tet_map=0;
  pyr_map=0;
  pri_map=0;
  hex_map=0;
  poly_map=0;
  if (num_procs > 1)
  {
    tri_map = (int**)malloc(nb*sizeof(int*));
    quad_map = (int**)malloc(nb*sizeof(int*));
    ngon_map = (int**)malloc(nb*sizeof(int*));
  }
  t_n = (int***)malloc(nb*sizeof(int**));
  q_n = (int***)malloc(nb*sizeof(int**));
  ngon_n = (int***)malloc(nb*sizeof(int**));
  for (b=0; b < nb; b++)
  {
    if (nt[b] > 0)
    {
      //fprintf(out_f,"\nBoundary %d, allocating space for %d triangles.",b,nt[b]);
      t_n[b] = (int**)malloc(nt[b]*sizeof(int*));
      for (i=0; i < nt[b]; i++)
        t_n[b][i] = (int*)malloc(3*sizeof(int));
      if (num_procs > 1)
        tri_map[b] = (int*)malloc(nt[b]*sizeof(int));
    } else
      t_n[b] = 0;
    if (nq[b] > 0)
    {
      //fprintf(out_f,"\nBoundary %d, allocating space for %d quadrilaterals.",b,nq[b]);
      q_n[b] = (int**)malloc(nq[b]*sizeof(int*));
      for (i=0; i < nq[b]; i++)
        q_n[b][i] = (int*)malloc(4*sizeof(int));
      if (num_procs > 1)
        quad_map[b] = (int*)malloc(nq[b]*sizeof(int));
    } else
      q_n[b] = 0;
    if (ng[b] > 0)
    {
      //fprintf(out_f,"\nBoundary %d, allocating space for %d polygons.",b,ng[b]);
      ngon_n[b] = (int**)malloc(ng[b]*sizeof(int*));
      for (i=0; i < ng[b]; i++)
        q_n[b][i] = 0;
      if (num_procs > 1)
        ngon_map[b] = (int*)malloc(ng[b]*sizeof(int));
    } else
      ngon_n[b] = 0;
  }
  if (ntet > 0)
  {
    //fprintf(out_f,"\nAllocating space for %d tetrahedra.",ntet);
    tet_n = (int**)malloc(ntet*sizeof(int*));
    for (i=0; i < ntet; i++)
      tet_n[i] = (int*)malloc(4*sizeof(int));
    if (num_procs > 1)
      tet_map = (int*)malloc(ntet*sizeof(int));
  } else
    tet_n = 0;
  if (npyr > 0)
  {
    //fprintf(out_f,"\nAllocating space for %d pyramids.",npyr);
    pyr_n = (int**)malloc(npyr*sizeof(int*));
    for (i=0; i < npyr; i++)
      pyr_n[i] = (int*)malloc(5*sizeof(int));
    if (num_procs > 1)
      pyr_map = (int*)malloc(npyr*sizeof(int));
  } else
    pyr_n = 0;
  if (npri > 0)
  {
    //fprintf(out_f,"\nAllocating space for %d prisms.",npri);
    pri_n = (int**)malloc(npri*sizeof(int*));
    for (i=0; i < npri; i++)
      pri_n[i] = (int*)malloc(6*sizeof(int));
    if (num_procs > 1)
      pri_map = (int*)malloc(npri*sizeof(int));
  } else
    pri_n = 0;
  if (nhex > 0)
  {
    //fprintf(out_f,"\nAllocating space for %d hexahedra.",nhex);
    hex_n = (int**)malloc(nhex*sizeof(int*));
    for (i=0; i < nhex; i++)
      hex_n[i] = (int*)malloc(8*sizeof(int));
    if (num_procs > 1)
      hex_map = (int*)malloc(nhex*sizeof(int));
  } else
    hex_n = 0;
  if (nply > 0)
  {
    //fprintf(out_f,"\nAllocating space for %d polyhedra.",nply);
    poly_n = (int***)malloc(nply*sizeof(int**));
    for (i=0; i < nply; i++)
      poly_n[i] = 0;
    if (num_procs > 1)
      poly_map = (int*)malloc(nply*sizeof(int));
  } else
    poly_n = 0;

  // store elements in temporary arrays
  ntet = npyr = npri = nhex = nply = 0;
  for (e=0; e < nelems; e++)
  {
    switch(element_type(e))
    {
      case 1:
        break;
      case 4:
        tet_n[ntet][0] = element[e].node_list->list[0];
        tet_n[ntet][1] = element[e].node_list->list[1];
        tet_n[ntet][2] = element[e].node_list->list[2];
        tet_n[ntet][3] = element[e].node_list->list[3];
        //f = element[e].face_list->list[0];
        //if (f > 0)
        //{
        //  tet_n[ntet][0] = face[f].node_list->list[0];
        //  tet_n[ntet][1] = face[f].node_list->list[2];
        //  tet_n[ntet][2] = face[f].node_list->list[1];
        //} else
        //{
        //  f = abs(f);
        //  tet_n[ntet][0] = face[f].node_list->list[0];
        //  tet_n[ntet][1] = face[f].node_list->list[1];
        //  tet_n[ntet][2] = face[f].node_list->list[2];
        //}
        //// find 4th node
        //m = abs(element[e].face_list->list[1]);
        //for (i=0; i < 3; i++)
        //  if (!face[f].node_list->Is_In_List(face[m].node_list->list[i]))
        //  {
        //    tet_n[ntet][3] = face[m].node_list->list[i];
        //    break;
        //  }
        ntet++;
        break;
      case 5:
        pyr_n[npyr][0] = element[e].node_list->list[0];
        pyr_n[npyr][1] = element[e].node_list->list[1];
        pyr_n[npyr][2] = element[e].node_list->list[2];
        pyr_n[npyr][3] = element[e].node_list->list[3];
        pyr_n[npyr][4] = element[e].node_list->list[4];
        //for (i=0; i < element[e].face_list->max; i++)
        //{
        //  f = element[e].face_list->list[i];
        //  if (face[abs(f)].node_list->max == 4)
        //  {
        //    if (f > 0)
        //    {
        //      pyr_n[npyr][0] = face[f].node_list->list[0];
        //      pyr_n[npyr][1] = face[f].node_list->list[3];
        //      pyr_n[npyr][2] = face[f].node_list->list[2];
        //      pyr_n[npyr][3] = face[f].node_list->list[1];
        //    } else
        //    {
        //      f = abs(f);
        //      pyr_n[npyr][0] = face[f].node_list->list[0];
        //      pyr_n[npyr][1] = face[f].node_list->list[1];
        //      pyr_n[npyr][2] = face[f].node_list->list[2];
        //      pyr_n[npyr][3] = face[f].node_list->list[3];
        //    }
        //    break;
        //  }
        //}
        //// find 5th node
        //pyr_n[npyr][4] = -1;
        //for (int l=0; l < element[e].face_list->max && pyr_n[npyr][4] < 0; l++)
        //{
        //  if (l==i) continue;
        //  int k = abs(element[e].face_list->list[l]);
        //  for (int m=0; m < face[k].node_list->max && pyr_n[npyr][4] < 0; m++)
        //    if (!face[f].node_list->Is_In_List(face[k].node_list->list[m]))
        //      pyr_n[npyr][4] = face[k].node_list->list[m];
        //}
        npyr++;
        break;
      case 6:
        pri_n[npri][0] = element[e].node_list->list[0];
        pri_n[npri][1] = element[e].node_list->list[1];
        pri_n[npri][2] = element[e].node_list->list[2];
        pri_n[npri][3] = element[e].node_list->list[3];
        pri_n[npri][4] = element[e].node_list->list[4];
        pri_n[npri][5] = element[e].node_list->list[5];
        //for (i=0; i < element[e].face_list->max; i++)
        //{
        //  f = element[e].face_list->list[i];
        //  if (face[abs(f)].node_list->max == 3)
        //  {
        //    if (f > 0)
        //    {
        //      pri_n[npri][0] = face[f].node_list->list[0];
        //      pri_n[npri][1] = face[f].node_list->list[2];
        //      pri_n[npri][2] = face[f].node_list->list[1];
        //    } else
        //    {
        //      f = abs(f);
        //      pri_n[npri][0] = face[f].node_list->list[0];
        //      pri_n[npri][1] = face[f].node_list->list[1];
        //      pri_n[npri][2] = face[f].node_list->list[2];
        //    }
        //    break;
        //  }
        //}
        //// now identify nodes above first three nodes in order
        //for (s=0; s < 3; s++)
        //{
        //  int m0, m1, m2;
        //  int n0 = pri_n[npri][s];
        //  int n1 = pri_n[npri][(s+1)%3];
        //  k=0;
        //  for (i=0; i < element[e].face_list->max && !k; i++)
        //  {
        //    f = element[e].face_list->list[i];
        //    if (face[abs(f)].node_list->max == 4)
        //    {
        //      for (j=0; j < face[abs(f)].node_list->max && !k; j++)
        //      {
        //        if (f > 0)
        //        {
        //          m0 = face[f].node_list->list[j];
        //          m1 = face[f].node_list->list[(j+1)%4];
        //          m2 = face[f].node_list->list[(j+3)%4];
        //        } else
        //        {
        //          m0 = face[-f].node_list->list[(j+1)%4];
        //          m1 = face[-f].node_list->list[j];
        //          m2 = face[-f].node_list->list[(j+2)%4];
        //        }
        //        if (m0 == n0 && m1 == n1)
        //        {
        //          pri_n[npri][s+3] = m2;
        //          k=1;
        //        }
        //      }
        //    }
        //  }
        //}
        npri++;
        break;
      case 8:
        hex_n[nhex][0] = element[e].node_list->list[0];
        hex_n[nhex][1] = element[e].node_list->list[1];
        hex_n[nhex][2] = element[e].node_list->list[2];
        hex_n[nhex][3] = element[e].node_list->list[3];
        hex_n[nhex][4] = element[e].node_list->list[4];
        hex_n[nhex][5] = element[e].node_list->list[5];
        hex_n[nhex][6] = element[e].node_list->list[6];
        hex_n[nhex][7] = element[e].node_list->list[7];
        //f = element[e].face_list->list[0];
        //if (f > 0)
        //{
        //  hex_n[nhex][0] = face[f].node_list->list[0];
        //  hex_n[nhex][1] = face[f].node_list->list[3];
        //  hex_n[nhex][2] = face[f].node_list->list[2];
        //  hex_n[nhex][3] = face[f].node_list->list[1];
        //} else
        //{
        //  f = abs(f);
        //  hex_n[nhex][0] = face[f].node_list->list[0];
        //  hex_n[nhex][1] = face[f].node_list->list[1];
        //  hex_n[nhex][2] = face[f].node_list->list[2];
        //  hex_n[nhex][3] = face[f].node_list->list[3];
        //}
        //// now process to get 4 remaining nodes in order
        //for (s=0; s < 4; s++)
        //{
        //  int m0, m1, m2;
        //  int n0 = hex_n[nhex][s];
        //  int n1 = hex_n[nhex][(s+1)%4];
        //  k=0;
        //  for (i=1; i < element[e].face_list->max && !k; i++)
        //  {
        //    f = element[e].face_list->list[i];
        //    for (j=0; j < face[abs(f)].node_list->max && !k; j++)
        //    {
        //      if (f > 0)
        //      {
        //        m0 = face[f].node_list->list[j];
        //        m1 = face[f].node_list->list[(j+1)%4];
        //        m2 = face[f].node_list->list[(j+3)%4];
        //      } else
        //      {
        //        m0 = face[-f].node_list->list[(j+1)%4];
        //        m1 = face[-f].node_list->list[j];
        //        m2 = face[-f].node_list->list[(j+2)%4];
        //      }
        //      if (m0 == n0 && m1 == n1)
        //      {
        //        hex_n[nhex][s+4] = m2;
        //        k=1;
        //      }
        //    }
        //  }
        //}
        nhex++;
        break;
      default:
        poly_n[nply] = (int**)malloc((element[e].face_list->max+1)*sizeof(int*));
        poly_n[nply][0] = (int*)malloc(sizeof(int));
        poly_n[nply][0][0] = element[e].face_list->max;
        for (i=0; i < element[e].face_list->max; i++)
        {
          f = element[e].face_list->list[i];
          poly_n[nply][i+1] = (int*)malloc((face[abs(f)].node_list->max+1)*sizeof(int));
          poly_n[nply][i+1][0] = face[abs(f)].node_list->max;
          if (f > 0)
          {
            for (j=0; j < face[f].node_list->max; j++)
              poly_n[nply][i+1][j+1] = face[f].node_list->list[j];
          } else
          {
            f = abs(f);
            for (k=1,j=face[f].node_list->max-1; j >= 0; j--)
              poly_n[nply][i+1][k++] = face[f].node_list->list[j];
          }
        }
        nply++;
        break;
    }
  }
  for (b=0; b < nb; b++)
  {
    nt[b] = nq[b] = ng[b] = 0;
    for (i=0; i < boundary[b].polygon_list->max; i++)
    {
      f = boundary[b].polygon_list->list[i];
      if (face[abs(f)].node_list->max == 3)
      {
        if (f > 0)
        {
          t_n[b][nt[b]][0] = face[f].node_list->list[0];
          t_n[b][nt[b]][1] = face[f].node_list->list[1];
          t_n[b][nt[b]][2] = face[f].node_list->list[2];
        } else
        {
          f = abs(f);
          t_n[b][nt[b]][0] = face[f].node_list->list[0];
          t_n[b][nt[b]][1] = face[f].node_list->list[2];
          t_n[b][nt[b]][2] = face[f].node_list->list[1];
        }
        nt[b]++;
      }
      if (face[abs(f)].node_list->max == 4)
      {
        if (f > 0)
        {
          q_n[b][nq[b]][0] = face[f].node_list->list[0];
          q_n[b][nq[b]][1] = face[f].node_list->list[1];
          q_n[b][nq[b]][2] = face[f].node_list->list[2];
          q_n[b][nq[b]][3] = face[f].node_list->list[3];
        } else
        {
          f = abs(f);
          q_n[b][nq[b]][0] = face[f].node_list->list[0];
          q_n[b][nq[b]][1] = face[f].node_list->list[3];
          q_n[b][nq[b]][2] = face[f].node_list->list[2];
          q_n[b][nq[b]][3] = face[f].node_list->list[1];
        }
        nq[b]++;
      }
      if (face[abs(f)].node_list->max > 4)
      {
        ngon_n[b][ng[b]] = (int*)malloc((face[abs(f)].node_list->max+1)*sizeof(int));
        ngon_n[b][ng[b]][0] = face[abs(f)].node_list->max;
        if (f > 0)
        {
          for (i=0; i < face[f].node_list->max; i++)
            ngon_n[b][ng[b]][i+1] = face[f].node_list->list[i];
        } else
        {
          f = abs(f);
          for (k=1,i=face[f].node_list->max-1; i >= 0; i--)
            ngon_n[b][ng[b]][k++] = face[f].node_list->list[i];
        }
        ng[b]++;
      }
    }
  }

  if (num_procs > 1)
  {
    create_hybrid_maps(nn,nb,ntet,npyr,npri,nhex,nply,tet_n,pyr_n,pri_n,hex_n,poly_n,
                       nt,t_n,nq,q_n,ng,ngon_n,nmap,tri_map,quad_map,ngon_map,
                       tet_map,pyr_map,pri_map,hex_map,poly_map);
    flag=P_CGNS_write(parallel,sname,nn,p,nb,bname,nt,t_n,nq,q_n,ng,ngon_n,
                      ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
                      nmap,tri_map,quad_map,ngon_map,tet_map,pyr_map,pri_map,hex_map,poly_map);
  } else
  {
    //no need for maps since no need to reassemble
    flag=CGNS_write(sname,nn,p,nb,bname,nt,t_n,nq,q_n,ng,ngon_n,ntet,tet_n,npyr,pyr_n,
                    npri,pri_n,nhex,hex_n,nply,poly_n); 
  }

  // clean memory
  free(p);
  if (num_procs > 1)
  {
    for (n=0; n < nn; n++)
      free(nmap[n]);
    free(nmap);
  }
  for (b=0; b < nb; b++)
  {
    free(bname[b]);
    if (nt[b] > 0)
    {
      for (i=0; i < nt[b]; i++)
        free(t_n[b][i]);
      free(t_n[b]);
      if (num_procs > 1)
        free(tri_map[b]);
    }
    if (nq[b] > 0)
    {
      for (i=0; i < nq[b]; i++)
        free(q_n[b][i]);
      free(q_n[b]);
      if (num_procs > 1)
        free(quad_map[b]);
    }
    if (ng[b] > 0)
    {
      for (i=0; i < ng[b]; i++)
        free(ngon_n[b][i]);
      free(ngon_n[b]);
      if (num_procs > 1)
        free(ngon_map[b]);
    }
  }
  free(bname);
  free(t_n);
  free(q_n);
  free(ngon_n);
  if (num_procs > 1)
  {
    free(tri_map);
    free(quad_map);
    free(ngon_map);
  }
  if (ntet > 0)
  {
    for (i=0; i < ntet; i++)
      free(tet_n[i]);
    free(tet_n);
    if (num_procs > 1)
      free(tet_map);
  }
  if (npyr > 0)
  {
    for (i=0; i < npyr; i++)
      free(pyr_n[i]);
    free(pyr_n);
    if (num_procs > 1)
      free(pyr_map);
  }
  if (npri > 0)
  {
    for (i=0; i < npri; i++)
      free(pri_n[i]);
    free(pri_n);
    if (num_procs > 1)
      free(pri_map);
  }
  if (nhex > 0)
  {
    for (i=0; i < nhex; i++)
      free(hex_n[i]);
    free(hex_n);
    if (num_procs > 1)
      free(hex_map);
  }
  if (nply > 0)
  {
    for (i=0; i < nply; i++)
    {
      for (j=poly_n[i][0][0]; j >= 0; j--)
        free(poly_n[i][j]);
      free(poly_n[i]);
    }
    free(poly_n);
    if (num_procs > 1)
      free(poly_map);
  }

  return (flag);
}

int POLYMESH::element_type(int e) // Identify standard element type
{
  int f, i, j, k, l, m, n;
  int type = 0;

  if (element[e].face_list->max == 1)
    type = 1;  // Shouldn't be used this way (flat elements are stored as polygons), but just in case!

  if (type == 0 && element[e].face_list->max == 4 && element[e].node_list->max == 4)
  {
    type = 4;
    for (i=0; i < element[e].face_list->max && type == 4; i++)
    {
      f = abs(element[e].face_list->list[i]);
      if (face[f].node_list->max != 3)
        type = 0;
    }
  }

  if (type == 0 && element[e].face_list->max == 5)
  {
    j = k = 0;
    for (i=0; i < element[e].face_list->max; i++)
    {
      f = abs(element[e].face_list->list[i]);
      if (face[f].node_list->max == 3) j++;
      if (face[f].node_list->max == 4) k++;
    }
    if (k==1 && j==4 && element[e].node_list->max == 5) type = 5;
    if (k==3 && j==2 && element[e].node_list->max == 6)
    {
      type = 6;
      // make sure triangles don't share nodes
      l = m = -1;
      for (i=0; i < element[e].face_list->max; i++)
      {
        f = abs(element[e].face_list->list[i]);
        if (face[f].node_list->max == 3)
          if (l < 0)
            l = f;
          else
            m = f;
      }
      for (j=0; j < face[l].node_list->max && type == 6; j++)
      {
        n = face[l].node_list->list[j];
        if (face[m].node_list->Is_In_List(n))
          type = 0;
      }
    }
  }

  if (type == 0 && element[e].face_list->max == 6 && element[e].node_list->max == 8)
  {
    type = 8;
    for (i=0; i < element[e].face_list->max && type == 8; i++)
    {
      f = abs(element[e].face_list->list[i]);
      if (face[f].node_list->max != 4)
        type = 0;
    }
  }

  return(type);
}

void POLYMESH::create_element_maps() // Identify owning proc and create global element number
{
  // THIS ASSUMES THE NODE OWNERSHIP IS CORRECT!!!
  int e, j, n, o;
#ifdef PARALLEL
  int b, f, i, k, l, m, p, ne, ln;
  MPI_Status single_status;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt;
  int *recvcnt;
  MPI_Request *srequest;
  MPI_Request *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff;
  char **rbuff;
  List **pelems;
#endif

  // first identify owning processors for each element as
  // lowest processor from any of the nodes of the element
  for (e=0; e < nelems; e++)
  {
    o = num_procs;
    for (j=0; j < element[e].node_list->max; j++)
    {
      n = element[e].node_list->list[j];
      o = MIN(o,node[n].me.proc);
    }
    //for (j=0; j < element[e].face_list->max; j++)
    //{
    //  f = abs(element[e].face_list->list[j]);
    //  for (i=0; i < face[f].node_list->max; i++)
    //  {
    //    n = face[f].node_list->list[i];
    //    o = MIN(o,node[n].me.proc);
    //  }
    //}
    if (o == my_rank)
      element[e].me = Pmap(-1,o,e);
    else
      element[e].me = Pmap(-1,o,-1);
  }

#ifdef PARALLEL
  // now take turns and set global element numbers
  if (num_procs > 1)
  {
    sendcnt = new int[num_procs];
    recvcnt = new int[num_procs];
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses  = new MPI_Status[num_procs];
    // allocate space for message buffers
    bdim = new int[num_procs];
    pelems = new List*[num_procs];
    for (p=0; p < num_procs; p++)
    {
      bdim[p] = 0;
      pelems[p] = new List();
    }
    sbuff = new char*[num_procs];
    rbuff = new char*[num_procs];
    for (p=0; p < num_procs; p++)
    {
      sbuff[p] = 0;
      rbuff[p] = 0;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int ntotal;
  // Assign global element numbers
  if (my_rank == 0)
  {
    ntotal = 0;
    for (e=0; e < nelems; e++)
      if (element[e].me.proc == my_rank)
        element[e].me.global = ntotal++;
    //fprintf(out_f,"\nCreate_Element_Maps: Proc %d, Current total number of elements = %d",my_rank,ntotal);
    //fflush(out_f);
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&ntotal,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&ntotal,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    fprintf(out_f,"\nCreate_Element_Maps: Total number of elements = %d",ntotal);
    fflush(out_f);
  } else
  {
#ifdef PARALLEL
    MPI_Recv(&ntotal,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < nelems; e++)
      if (element[e].me.proc == my_rank)
        element[e].me.global = ntotal++;
    MPI_Send(&ntotal,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
#endif
    //fprintf(out_f,"\nCreate_Element_Maps: Proc %d, Current total number of elements = %d",my_rank,ntotal);
    //fflush(out_f);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);

  if (num_procs > 1)
  {
    // communicate and exchange the element global number and indices not owned by me
    List **elist;

    elist = new List*[nn];

    for (n=0; n < nn; n++)
      elist[n] = new List(); // hash table (list of elements per node)

    for (e=0; e < nelems; e++)
    {
      for (j=0; j < element[e].node_list->max; j++)
      {
        n = element[e].node_list->list[j];
        elist[n]->Add_To_List(e);
      }
      //for (j=0; j < element[e].face_list->max; j++)
      //{
      //  f = abs(element[e].face_list->list[j]);
      //  for (i=0; i < face[f].node_list->max; i++)
      //  {
      //    n = face[f].node_list->list[i];
      //    elist[n]->Check_List(e);
      //  }
      //}
    }

    // create list of elements to send to each processor
    for (e=0; e < nelems; e++)
      if ((p=element[e].me.proc) != my_rank)
        pelems[p]->Add_To_List(e);

    // send buffer size
    for (p=0; p < num_procs; p++)
    {
      sendcnt[p] = 0;
      if (p != my_rank)
      {
        sendcnt[p] = sizeof(int);
        for (i=0; i < pelems[p]->max; i++)
        {
          e = pelems[p]->list[i];
          sendcnt[p] += (element[e].node_list->max*3 + 2)*sizeof(int);
        }
      }
    }
  
    nreq_s = nreq_r = 0;
    for (p=0; p < num_procs; p++)
    {
      recvcnt[p] = 0;
      if (p != my_rank)
      {
        MPI_Isend(&sendcnt[p],1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_s++;
        nreq_r++;
      }
    }
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    // size message buffers for each processor
    for (p=0; p < num_procs; p++)
    {
      bdim[p] = 0;
      if (p == my_rank || (sendcnt[p] == 0 && recvcnt[p] == 0))
        continue;

      // size both buffers to p_i for voxel nodes
      bdim[p] = MAX(sendcnt[p],recvcnt[p]);
      if (bdim[p] > 0) sbuff[p] = new char[bdim[p]];
      if (bdim[p] > 0) rbuff[p] = new char[bdim[p]];
    }

    // package the list of element per processor
    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank || sendcnt[p] == 0)
        continue;

      sposition=0;
      MPI_Pack(&(pelems[p]->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      for (i=0; i < pelems[p]->max; i++)
      {
        e = pelems[p]->list[i];
        MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(element[e].node_list->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        for (j=0; j < element[e].node_list->max; j++)
        {
          n = element[e].node_list->list[j];
          MPI_Pack(&(node[n].me.global),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          MPI_Pack(&(node[n].me.proc),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          MPI_Pack(&(node[n].me.index),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        }
      }
    }

    // exchange messages
    nreq_s = nreq_r = 0;
    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;
      if (sendcnt[p] > 0)
      {
        MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }
      if (recvcnt[p] > 0)
      {
        MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    // unpack element list and re-package the current map for requested elements
    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;

      sposition=rposition=0;
      Pmap *nmap;
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Pack(&ne,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      for (i=0; i < ne; i++)
      {
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        int esave = e;
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&ln,1,MPI_INT,MPI_COMM_WORLD);
        nmap = new Pmap[ln];
        m = -1;
        for (n=0; n < ln; n++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&nmap[n].global,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&nmap[n].proc,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&nmap[n].index,1,MPI_INT,MPI_COMM_WORLD);
          // find a node to use with the local node-element hash table
          if (m < 0 && nmap[n].proc == my_rank)
            m = nmap[n].index;
        }
        if (m < 0)
        {
          fprintf(out_f,"\nCreate_Element_Maps: no local node identified!");
          fflush(out_f);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        }

        int le = -1; // local element identity
        for (j=0; j < elist[m]->max && le < 0; j++)
        {
          e = elist[m]->list[j];
          if (element[e].node_list->max != ln)
            continue;
          bool ematch = true;
          for (k=0; k < ln && ematch; k++)
          {
            n = element[e].node_list->list[k];
            bool nmatch = false;
            for (l=0; l < ln && !nmatch; l++)
              if (nmap[l] == node[n].me)
                nmatch = true;
            ematch = nmatch;
          }
          if (ematch)
            le = e;
        }

        if (le < 0)
        {
          fprintf(out_f,"\nCreate_Element_Maps: proc %d no local element identified!",my_rank);
          fprintf(out_f,"\n Proc %d, i = %d, passed element = %d from process %d",my_rank,i,esave,p);
          for (n=0; n < ln; n++)
            fprintf(out_f,"\nProc %d, node %d = nmap = (%d, %d, %d)",my_rank,n,nmap[n].global,nmap[n].proc,nmap[n].index);
          fprintf(out_f,"\nProc %d, local node = %d",my_rank,m);
          for (j=0; j < elist[m]->max; j++)
          {
            e = elist[m]->list[j];
            fprintf(out_f,"\nProc %d, element %d = %d",my_rank,j,e);
            for (k=0; k < element[e].node_list->max; k++)
            {
              n = element[e].node_list->list[k];
              fprintf(out_f,"\nProc %d, element node %d = (%d, %d, %d)",my_rank,k,node[n].me.global,node[n].me.proc,node[n].me.index);
            }
          }
          fflush(out_f);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        }

        if (element[le].me.proc != my_rank)
        {
          fprintf(out_f,"\nCreate_Element_Maps: local element not owned by proc!");
          fflush(out_f);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        }
        
        MPI_Pack(&element[le].me.global,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&element[le].me.proc,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&element[le].me.index,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);

        delete[] nmap;
      }
    }

    // exchange messages
    nreq_s = nreq_r = 0;
    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;
      if (recvcnt[p] > 0)
      {
        MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }
      if (sendcnt[p] > 0)
      {
        MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    // unpack element list and re-package the current map for requested elements
    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;

      sposition=rposition=0;
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
      for (m=0; m < ne; m++)
      {
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&j,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&k,1,MPI_INT,MPI_COMM_WORLD);
        element[e].me = Pmap(i,j,k);
      }
    }

    for (n=0; n < nn; n++)
      delete elist[n];
    delete[] elist;

    // release communication memory
    for (p=0; p < num_procs; p++)
    {
      delete pelems[p];
      if (sbuff[p] > 0) delete[] sbuff[p];
      if (rbuff[p] > 0) delete[] rbuff[p];
    }
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] srequest;
    delete[] rrequest;
    delete[] statuses;
    // allocate space for message buffers
    delete[] bdim;
    delete[] pelems;
    delete[] sbuff;
    delete[] rbuff;
  }
#endif

  return;
}

