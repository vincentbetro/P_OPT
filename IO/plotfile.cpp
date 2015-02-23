#ifdef PARALLEL
#include "mpi.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Point.h"
#include "Vector.h"
#include "smooth.h"

extern FILE *in_f, *jou_f, *out_f; /*input output journal files global*/

extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/

//void POLYMESH::Ensight_layers(char sname[], int ne, int **edge)
//{
//  FILE *fp, *casefp, **funcp;
//  int const bdim = 80;
//  char buff[bdim], gname[bdim], (*fname)[bdim];
//  char filename[bdim], buff2[bdim], write_mode[4];


//  return;
//}

void POLYMESH::Ensight(int mode, char sname[], int num_func, char (*func_name)[80], double **func, int nl, int *ne, int ***edge)
{
  FILE *fp, *casefp, **funcp;
  int const bdim = 80;
  char buff[bdim], gname[bdim], (*fname)[bdim];
  char filename[bdim], buff2[bdim], write_mode[4];
  int b, e, f, i, j, k, l, n, part, proc;
  int ntet, npyr, npri, nhex, nply, nt, nq;
  float fx, fy, fz, ff;
  bool ascii = false;

  if (mode != 1) ascii = true;

  sprintf(buff,"%s",sname);
  if (num_procs > 1)
    sprintf(buff2,"_%d.cgns",my_rank);
  else
    sprintf(buff2,".cgns");
  char *ptr = strstr(buff,buff2);
  if (ptr == NULL)
  {
    fprintf(out_f,"\nCGNS suffix <.cgns> not found in file name!");
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
  if (nl > 0)
    strcat(buff,".layers");

  //
  // /write Ensight files
  //
  sprintf(filename,"%s.case",buff);
  if (my_rank == 0)
  {
    // Open file for write
    if ((casefp = fopen(filename,"w")) == 0)
    {
      printf("\nEnsight: Error opening file <%s>.",filename);
      exit(0);
    }

    fprintf(casefp,"# Ensight Gold Case file: %s\n\n",filename);

    fprintf(casefp,"FORMAT\n\n");
    fprintf(casefp,"type:  ensight gold\n\n");
  }

  sprintf(buff,"%s",filename);
  ptr = strstr(buff,".case");
  if (ptr == NULL)
  {
    printf("\nEnsight suffix <.case> not found in file name!");
    fflush(stdout);
    exit(0);
  } else
    *ptr = '\0';

  sprintf(gname,"%s.geo",buff);
  if (num_func > 0)
  {
    fname = new char[num_func][bdim];
    for (l=0; l < num_func; l++)
      sprintf(fname[l],"%s.%s",buff,func_name[l]);
    funcp = new FILE*[num_func];
  }

  if (my_rank == 0)
  {
    fprintf(casefp,"GEOMETRY\n\n");
    fprintf(casefp,"model:     %s\n\n",gname);

    fprintf(casefp,"VARIABLE\n\n");
    for (l=0; l < num_func; l++)
      fprintf(casefp,"scalar per node: %s %s\n\n",func_name[l],fname[l]);

    fclose(casefp);
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  part = 0;

  for (proc = 0; proc < num_procs; proc++)
  {

#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    i = part;
    MPI_Allreduce(&i,&part,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);        
#endif

    if (proc != my_rank)
      continue;

    // Open file for write
    if (ascii)
    {
      if (proc == 0)
        sprintf(write_mode,"w");
      else
        sprintf(write_mode,"a");
      if ((fp = fopen(gname,write_mode)) == 0)
      {
        printf("\nEnsight: Error opening geometry file <%s>.",gname);
        exit(0);
      }
      for (l=0; l < num_func; l++)
      {
        if ((funcp[l] = fopen(fname[l],write_mode)) == 0)
        {
          printf("\nEnsight: Error opening function file <%s>.",fname[l]);
          exit(0);
        }
      }
    } else
    {
      if (proc == 0)
        sprintf(write_mode,"wb");
      else
        sprintf(write_mode,"ab");
      if ((fp = fopen(gname,write_mode)) == 0)
      {
        printf("\nEnsight: Error opening geometry file <%s>.",gname);
        exit(0);
      }
      for (l=0; l < num_func; l++)
      {
        if ((funcp[l] = fopen(fname[l],write_mode)) == 0)
        {
          printf("\nEnsight: Error opening function file <%s>.",fname[l]);
          exit(0);
        }
      }
    }

    if (proc == 0)
    {
      buff[sizeof(buff)-1] = '\0';
      if (ascii)
      {
        fprintf(fp,"Ensight Geometry file written from P_SMOOTH.\n");
        fprintf(fp,"Case file name = %s\n",filename);
        for (l=0; l < num_func; l++)
          fprintf(funcp[l],"Function file %s written from P_SMOOTH.\n",func_name[l]);
      } else
      {
        strncpy(buff,"C Binary",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        strncpy(buff,"Ensight Geometry file written from P_SMOOTH.",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        sprintf(buff,"Case file name = %s",filename);
        fwrite(buff,sizeof(char),80,fp);
        for (l=0; l < num_func; l++)
        {
          sprintf(buff2,"Function file %s written from P_SMOOTH",func_name[l]);
          strncpy(buff,buff2,sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,funcp[l]);
        }
      }

      if (ascii)
      {
        fprintf(fp,"node id given\n");
        fprintf(fp,"element id given\n");
      } else
      {
        strncpy(buff,"node id given",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        strncpy(buff,"element id given",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
      }
    }

    if (ascii)
    {
      fprintf(fp,"part\n");
      for (l=0; l < num_func; l++)
        fprintf(funcp[l],"part\n");
      ++part;
      fprintf(fp,"%10d\n",part);
      for (l=0; l < num_func; l++)
        fprintf(funcp[l],"%10d\n",part);
      if (num_procs > 1)
        fprintf(fp,"Proc %d 3D uns-elements\n",proc);
      else
        fprintf(fp,"3D uns-elements\n");
    } else
    {
      strncpy(buff,"part",sizeof(buff)-1);
      fwrite(buff,sizeof(char),80,fp);
      for (l=0; l < num_func; l++)
        fwrite(buff,sizeof(char),80,funcp[l]);
      ++part;
      fwrite(&part,sizeof(int),1,fp);
      for (l=0; l < num_func; l++)
        fwrite(&part,sizeof(int),1,funcp[l]);
      if (num_procs > 1)
        sprintf(buff2,"Proc %d 3D uns-elements",proc);
      else
        sprintf(buff2,"3D uns-elements");
      strncpy(buff,buff2,sizeof(buff)-1);
      fwrite(buff,sizeof(char),80,fp);
    }

    int *map = new int[nn];
    for (n=0; n < nn; n++)
      map[n] = 0;
    for (e=0; e < nelems; e++)
    {
      if (element[e].me.proc != my_rank)
        continue;
      for (i=0; i < element[e].node_list->max; i++)
        map[element[e].node_list->list[i]] = 1;
    }
    i=0;
    for (n=0; n < nn; n++)
      if (map[n] > 0) map[n] = ++i;

    if (ascii)
    {
      fprintf(fp,"coordinates\n");
      i=0;
      for (n=0; n < nn; n++)
        if (map[n] > 0) i++;
      fprintf(fp,"%10d\n",i);
      for (n=0; n < nn; n++)
        if (map[n] > 0)
          fprintf(fp,"%10d\n",map[n]);
      for (n=0; n < nn; n++)
        if (map[n] > 0)
          fprintf(fp,"%12.5e\n",node[n].vert[0]);
      for (n=0; n < nn; n++)
        if (map[n] > 0)
          fprintf(fp,"%12.5e\n",node[n].vert[1]);
      for (n=0; n < nn; n++)
        if (map[n] > 0)
          fprintf(fp,"%12.5e\n",node[n].vert[2]);
      for (l=0; l < num_func; l++)
      {
        fprintf(funcp[l],"coordinates\n");
        for (n=0; n < nn; n++)
          if (map[n] > 0)
            fprintf(funcp[l],"%12.5e\n",func[l][n]);
      }
    } else
    {
      strncpy(buff,"coordinates",sizeof(buff)-1);
      fwrite(buff,sizeof(char),80,fp);
      for (l=0; l < num_func; l++)
        fwrite(buff,sizeof(char),80,funcp[l]);
      i=0;
      for (n=0; n < nn; n++)
        if (map[n] > 0) i++;
      fwrite(&i,sizeof(int),1,fp);
      for (n=0; n < nn; n++)
        if (map[n] > 0)
          fwrite(&map[n],sizeof(int),1,fp);
      for (n=0; n < nn; n++)
      {
        if (map[n] > 0)
        {
          fx = (float)node[n].vert[0];
          fwrite(&fx,sizeof(float),1,fp);
        }
      }
      for (n=0; n < nn; n++)
      {
        if (map[n] > 0)
        {
          fy = (float)node[n].vert[1];
          fwrite(&fy,sizeof(float),1,fp);
        }
      }
      for (n=0; n < nn; n++)
      {
        if (map[n] > 0)
        {
          fz = (float)node[n].vert[2];
          fwrite(&fz,sizeof(float),1,fp);
        }
      }
      if (num_func > 0)
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            for (l=0; l < num_func; l++)
            {
              ff = (float)func[l][n];
              fwrite(&ff,sizeof(float),1,funcp[l]);
            }
          }
        }
    }

    ntet = npyr = npri = nhex = nply = 0;
    for (e=0; e < nelems; e++)
    {
      if (element[e].me.proc != my_rank)
        continue;
      switch(element_type(e))
      {
        case 4: ntet++; break;
        case 5: npyr++; break;
        case 6: npri++; break;
        case 8: nhex++; break;
        default: nply++; break;
      }
    }

    if (ntet > 0)
    {
      if (ascii)
      {
        fprintf(fp,"tetra4\n");
        fprintf(fp,"%10d\n",ntet);
        for (n=0; n < ntet; n++)
          fprintf(fp,"%10d\n",n+1);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 4) continue;
          for (j=0; j < 4; j++)
            fprintf(fp,"%10d",map[element[e].node_list->list[j]]);
          fprintf(fp,"\n");
        }
      } else
      {
        strncpy(buff,"tetra4",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&ntet,sizeof(int),1,fp);
        for (n=1; n <= ntet; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 4) continue;
          for (j=0; j < 4; j++)
          {
            i=map[element[e].node_list->list[j]];
            fwrite(&i,sizeof(int),1,fp);
          }
        }
      }
    }
    if (npyr > 0)
    {
      if (ascii)
      {
        fprintf(fp,"pyramid5\n");
        fprintf(fp,"%10d\n",npyr);
        for (n=0; n < npyr; n++)
          fprintf(fp,"%10d\n",n+1);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 5) continue;
          for (j=0; j < 5; j++)
            fprintf(fp,"%10d",map[element[e].node_list->list[j]]);
          fprintf(fp,"\n");
        }
      } else
      {
        strncpy(buff,"pyramid5",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&npyr,sizeof(int),1,fp);
        for (n=1; n <= npyr; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 5) continue;
          for (j=0; j < 5; j++)
          {
            i=map[element[e].node_list->list[j]];
            fwrite(&i,sizeof(int),1,fp);
          }
        }
      }
    }
    if (npri > 0)
    {
      if (ascii)
      {
        fprintf(fp,"penta6\n");
        fprintf(fp,"%10d\n",npri);
        for (n=0; n < npri; n++)
          fprintf(fp,"%10d\n",n+1);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 6) continue;
          for (j=0; j < 6; j++)
            fprintf(fp,"%10d",map[element[e].node_list->list[j]]);
          fprintf(fp,"\n");
        }
      } else
      {
        strncpy(buff,"penta6",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&npri,sizeof(int),1,fp);
        for (n=1; n <= npri; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 6) continue;
          for (j=0; j < 6; j++)
          {
            i=map[element[e].node_list->list[j]];
            fwrite(&i,sizeof(int),1,fp);
          }
        }
      }
    }
    if (nhex > 0)
    {
      if (ascii)
      {
        fprintf(fp,"hexa8\n");
        fprintf(fp,"%10d\n",nhex);
        for (n=0; n < nhex; n++)
          fprintf(fp,"%10d\n",n+1);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 8) continue;
          for (j=0; j < 8; j++)
            fprintf(fp,"%10d",map[element[e].node_list->list[j]]);
          fprintf(fp,"\n");
        }
      } else
      {
        strncpy(buff,"hexa8",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&nhex,sizeof(int),1,fp);
        for (n=1; n <= nhex; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          if (element_type(e) != 8) continue;
          for (j=0; j < 8; j++)
          {
            i=map[element[e].node_list->list[j]];
            fwrite(&i,sizeof(int),1,fp);
          }
        }
      }
    }
    if (nply > 0)
    {
      if (ascii)
      {
        fprintf(fp,"nfaced\n");
        fprintf(fp,"%10d\n",nply);
        for (n=0; n < nply; n++)
          fprintf(fp,"%10d\n",n+1);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          fprintf(fp,"%10d\n",element[e].face_list->max);
        }
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          for (j=0; j < element[e].face_list->max; j++)
          {
            f = abs(element[e].face_list->list[j]);
            fprintf(fp,"%10d\n",face[f].node_list->max);
          }
        }
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          for (j=0; j < element[e].face_list->max; j++)
          {
            f = element[e].face_list->list[j];
            if (f < 0)
            {
              for (k=face[abs(f)].node_list->max-1; k >= 0; k--)
                fprintf(fp,"%10d",map[face[abs(f)].node_list->list[k]]);
              fprintf(fp,"\n");
            } else
            {
              for (k=0; k < face[f].node_list->max; k++)
                fprintf(fp,"%10d",map[face[f].node_list->list[k]]);
              fprintf(fp,"\n");
            }
          }
        }
      } else
      {
        strncpy(buff,"nfaced",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&nply,sizeof(int),1,fp);
        for (n=1; n <= nply; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          fwrite(&(element[e].face_list->max),sizeof(int),1,fp);
        }
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          for (j=0; j < element[e].face_list->max; j++)
          {
            f = abs(element[e].face_list->list[j]);
            fwrite(&(face[f].node_list->max),sizeof(int),1,fp);
          }
        }
        for (e=0; e < nelems; e++)
        {
          if (element[e].me.proc != my_rank)
            continue;
          k= element_type(e);
          if (k==4 || k == 5 || k == 6 || k == 8) continue;
          for (j=0; j < element[e].face_list->max; j++)
          {
            f = element[e].face_list->list[j];
            if (f < 0)
            {
              for (k=face[abs(f)].node_list->max-1; k >= 0; k--)
                fwrite(&map[face[abs(f)].node_list->list[k]],sizeof(int),1,fp);
            } else
            {
              for (k=0; k < face[f].node_list->max; k++)
                fwrite(&map[face[f].node_list->list[k]],sizeof(int),1,fp);
            }
          }
        }
      }
    }
    
    for (b=0; b < nb; b++)
    {
      for (n=0; n < nn; n++)
        map[n] = 0;
      for (i=0; i < boundary[b].polygon_list->max; i++)
      {
        f = abs(boundary[b].polygon_list->list[i]);
        for (j=0; j < face[f].node_list->max; j++)
          map[face[f].node_list->list[j]] = 1;
      }
      i=0;
      for (n=0; n < nn; n++)
        if (map[n] > 0) map[n] = ++i;

      if (i == 0)
        continue;

      if (ascii)
      {
        ++part;
        fprintf(fp,"part\n%10d\n",part);
        for (l=0; l < num_func; l++)
          fprintf(funcp[l],"part\n%10d\n",part);
        if (num_procs > 1)
          fprintf(fp,"Proc %d %s\n",proc,boundary[b].name);
        else
          fprintf(fp,"%s\n",boundary[b].name);
      } else
      {
        strncpy(buff,"part",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        for (l=0; l < num_func; l++)
          fwrite(buff,sizeof(char),80,funcp[l]);
        ++part;
        fwrite(&part,sizeof(int),1,fp);
        for (l=0; l < num_func; l++)
          fwrite(&part,sizeof(int),1,funcp[l]);
        strncpy(buff," ",sizeof(buff)-1);
        if (num_procs > 1)
          sprintf(buff,"Proc %d %s",proc,boundary[b].name);
        else
          sprintf(buff,"%s",boundary[b].name);
        fwrite(buff,sizeof(char),80,fp);
      }

      if (ascii)
      {
        fprintf(fp,"coordinates\n");
        fprintf(fp,"%10d\n",i);
        for (n=0; n < nn; n++)
          if (map[n] > 0)
            fprintf(fp,"%10d\n",map[n]);
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[0]);
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[1]);
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[2]);
        }
        for (l=0; l < num_func; l++)
        {
          fprintf(funcp[l],"coordinates\n");
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fprintf(funcp[l],"%12.5e\n",func[l][n]);
        }
      } else
      {
        strncpy(buff,"coordinates",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        for (l=0; l < num_func; l++)
          fwrite(buff,sizeof(char),80,funcp[l]);
        fwrite(&i,sizeof(int),1,fp);
        for (n=0; n < nn; n++)
          if (map[n] > 0)
            fwrite(&map[n],sizeof(int),1,fp);
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fx = (float)node[n].vert[0];
            fwrite(&fx,sizeof(float),1,fp);
          }
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fy = (float)node[n].vert[1];
            fwrite(&fy,sizeof(float),1,fp);
          }
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fz = (float)node[n].vert[2];
            fwrite(&fz,sizeof(float),1,fp);
          }
        }
        if (num_func > 0)
          for (n=0; n < nn; n++)
          {
            if (map[n] > 0)
            {
              for (l=0; l < num_func; l++)
              {
                ff = (float)func[l][n];
                fwrite(&ff,sizeof(float),1,funcp[l]);
              }
            }
          }
      }

      int ngon = 0;
      nt = nq = 0;
      for (i=0; i < boundary[b].polygon_list->max; i++)
      {
        f = abs(boundary[b].polygon_list->list[i]);

        if (face[f].node_list->max == 3)
          nt++;
        else if (face[f].node_list->max == 4)
          nq++;
        else
          ngon++;
      }

      if (nt > 0)
      {
        if (ascii)
        {
          fprintf(fp,"tria3\n");
          fprintf(fp,"%10d\n",nt);
          for (n=0; n < nt; n++)
            fprintf(fp,"%10d\n",n+1);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max != 3) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fprintf(fp,"%10d",map[face[abs(f)].node_list->list[j]]);
              fprintf(fp,"\n");
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fprintf(fp,"%10d",map[face[f].node_list->list[j]]);
              fprintf(fp,"\n");
            }
          }
        } else
        {
          strncpy(buff,"tria3",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&nt,sizeof(int),1,fp);
          for (n=1; n <= nt; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max != 3) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fwrite(&map[face[abs(f)].node_list->list[j]],sizeof(int),1,fp);
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fwrite(&map[face[f].node_list->list[j]],sizeof(int),1,fp);
            }
          }
        }
      }
      if (nq > 0)
      {
        if (ascii)
        {
          fprintf(fp,"quad4\n");
          fprintf(fp,"%10d\n",nq);
          for (n=0; n < nq; n++)
            fprintf(fp,"%10d\n",n+1);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max != 4) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fprintf(fp,"%10d",map[face[abs(f)].node_list->list[j]]);
              fprintf(fp,"\n");
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fprintf(fp,"%10d",map[face[f].node_list->list[j]]);
              fprintf(fp,"\n");
            }
          }
        } else
        {
          strncpy(buff,"quad4",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&(nq),sizeof(int),1,fp);
          for (n=1; n <= nq; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max != 4) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fwrite(&map[face[abs(f)].node_list->list[j]],sizeof(int),1,fp);
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fwrite(&map[face[f].node_list->list[j]],sizeof(int),1,fp);
            }
          }
        }
      }
      if (ngon > 0)
      {
        if (ascii)
        {
          fprintf(fp,"nsided\n");
          fprintf(fp,"%10d\n",ngon);
          for (n=0; n < ngon; n++)
            fprintf(fp,"%10d\n",n+1);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = abs(boundary[b].polygon_list->list[i]);
            if (face[f].node_list->max <= 4) continue;
            fprintf(fp,"%10d\n",face[f].node_list->max);
          }
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max <= 4) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fprintf(fp,"%10d",map[face[abs(f)].node_list->list[j]]);
              fprintf(fp,"\n");
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fprintf(fp,"%10d",map[face[f].node_list->list[j]]);
              fprintf(fp,"\n");
            }
          }
        } else
        {
          strncpy(buff,"nsided",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&(ngon),sizeof(int),1,fp);
          for (n=1; n <= ngon; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = abs(boundary[b].polygon_list->list[i]);
            if (face[f].node_list->max <= 4) continue;
            fwrite(&(face[f].node_list->max),sizeof(int),1,fp);
          }
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            if (face[abs(f)].node_list->max <= 4) continue;
            if (f < 0)
            {
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                fwrite(&map[face[abs(f)].node_list->list[j]],sizeof(int),1,fp);
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fwrite(&map[face[f].node_list->list[j]],sizeof(int),1,fp);
            }
          }
        }
      }
    }
    for (l=0; l < nl; l++)
    {
      for (n=0; n < nn; n++)
        map[n] = 0;
      for (i=0; i < ne[l]; i++)
      {
        for (j=0; j < 2; j++)
          map[edge[l][i][j]] = 1;
      }
      i=0;
      for (n=0; n < nn; n++)
        if (map[n] > 0) map[n] = ++i;

      if (i == 0)
        continue;

      if (ascii)
      {
        ++part;
        fprintf(fp,"part\n%10d\n",part);
        if (num_procs > 1)
          fprintf(fp,"Proc %d Layer %i\n",proc,l+1);
        else
          fprintf(fp,"Layer %i\n",l+1);
      } else
      {
        strncpy(buff,"part",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        ++part;
        fwrite(&part,sizeof(int),1,fp);
        strncpy(buff," ",sizeof(buff)-1);
        if (num_procs > 1)
          sprintf(buff,"Proc %d Layer %i",proc,l+1);
        else
          sprintf(buff,"Layer %i",l+1);
        fwrite(buff,sizeof(char),80,fp);
      }

      if (ascii)
      {
        fprintf(fp,"coordinates\n");
        fprintf(fp,"%10d\n",i);
        for (n=0; n < nn; n++)
          if (map[n] > 0)
            fprintf(fp,"%10d\n",map[n]);
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[0]);
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[1]);
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
            fprintf(fp,"%12.5e\n",node[n].vert[2]);
        }
      } else
      {
        strncpy(buff,"coordinates",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&i,sizeof(int),1,fp);
        for (n=0; n < nn; n++)
          if (map[n] > 0)
            fwrite(&map[n],sizeof(int),1,fp);
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fx = (float)node[n].vert[0];
            fwrite(&fx,sizeof(float),1,fp);
          }
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fy = (float)node[n].vert[1];
            fwrite(&fy,sizeof(float),1,fp);
          }
        }
        for (n=0; n < nn; n++)
        {
          if (map[n] > 0)
          {
            fz = (float)node[n].vert[2];
            fwrite(&fz,sizeof(float),1,fp);
          }
        }
      }

      if (ne[l] > 0)
      {
        if (ascii)
        {
          fprintf(fp,"bar2\n");
          fprintf(fp,"%10d\n",ne[l]);
          for (n=0; n < ne[l]; n++)
            fprintf(fp,"%10d\n",n+1);
          for (i=0; i < ne[l]; i++)
          {
            for (j=0; j < 2; j++)
              fprintf(fp,"%10d",map[edge[l][i][j]]);
            fprintf(fp,"\n");
          }
        } else
        {
          strncpy(buff,"bar2",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&ne[l],sizeof(int),1,fp);
          for (n=1; n <= ne[l]; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (i=0; i < ne[l]; i++)
          {
            for (j=0; j < 2; j++)
              fwrite(&map[edge[l][i][j]],sizeof(int),1,fp);
          }
        }
      }
    }

    delete[] map;

    fclose(fp);
    for (l=0; l < num_func; l++)
      fclose(funcp[l]);

  }

  if (num_func > 0)
  {
    delete[] fname;
    delete[] funcp;
  }

  return;

}

// Fieldview routines
#include "FV_routines.h"

void POLYMESH::Fieldview(int mode, char sname[], int num_func, char (*func_name)[80], double **func)
{
  FILE* fp;
  int const bdim = 80;
  char buff[bdim], buff2[bdim];
  char filename[bdim];
  int b, e, f, i, j, k, n;
  int n0, n1, n2, n3, n4, n5, n6, n7, proc;
  int csize, isize, fsize;
  int ntet, npyr, npri, nhex, nply, nt, nq;
  int ibuf[257];
  float *ftmp;

  csize = sizeof(char);
  isize = sizeof(int);
  fsize = sizeof(float);

  sprintf(buff,"%s",sname);
  if (num_procs > 1)
    sprintf(buff2,"_%d.cgns",my_rank);
  else
    sprintf(buff2,".cgns");
  char *ptr = strstr(buff,buff2);
  if (ptr == NULL)
  {
    fprintf(out_f,"\nCGNS suffix <.cgns> not found in file name!");
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

  //
  // write Fieldview file
  //
  switch(mode)
  {
    case 1: // Write Binary
      sprintf(filename,"%s.fv",buff);
      for (proc = 0; proc < num_procs; proc++)
      {
        #ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        if (proc != my_rank)
          continue;

        // Open file for write
        if (proc == 0)
        {
          if ((fp = fopen(filename,"wb")) == 0)
          {
            fprintf(out_f,"\nError opening file <%s>.",filename);
            #ifdef PARALLEL
            MPI_Abort(MPI_COMM_WORLD,0);
            #endif
            exit(0);
          }

          // Fieldview bit pattern
          ibuf[0] = FV_MAGIC;
          fwrite(ibuf,isize,1,fp);

          // name and version number
          sprintf(buff,"FIELDVIEW");
          fwrite(buff,csize,80,fp);
          ibuf[0] = 3;
          ibuf[1] = 0;
          fwrite(ibuf,isize,2,fp);
          ibuf[0] = FV_COMBINED_FILE;
          fwrite(ibuf,isize,1,fp);

          ibuf[0] = 0;
          fwrite(ibuf,isize,1,fp);

          ftmp = (float*)malloc(4*fsize);
          ftmp[0] = 1.0;
          ftmp[1] = 0.0;
          ftmp[2] = 0.0;
          ftmp[3] = 0.0;
          fwrite(ftmp,fsize,4,fp);
          free(ftmp);

          // Number of grids
          ibuf[0] = num_procs;
          fwrite(ibuf,isize,1,fp);

          // Number of boundary types
          ibuf[0] = nb;
          fwrite(ibuf,isize,1,fp);

          // Boundary types
          for (b=0; b < nb; b++)
          {
            for (i=0; i < 80; i++)
              buff[i] = ' ';
            sprintf(buff,"%s",boundary[b].name);
            ibuf[0] = 0;
            ibuf[1] = 1;
            fwrite(ibuf,isize,2,fp);
            fwrite(buff,csize,80,fp);
          }

          // Variable names
          ibuf[0] = num_func;
          fwrite(ibuf,isize,1,fp);
          for (n=0; n < num_func; n++)
          {
            for (i=0; i < 80; i++)
              buff[i] = ' ';
            sprintf(buff,"%s",func_name[n]);
            fwrite(buff,csize,80,fp);
          }

          // Boundary Variable names
          ibuf[0] = 0;
          fwrite(ibuf,isize,1,fp);
        } else
        {
          if ((fp = fopen(filename,"ab")) == 0)
          {
            fprintf(out_f,"\nError opening file <%s>.",filename);
            #ifdef PARALLEL
            MPI_Abort(MPI_COMM_WORLD,0);
            #endif
            exit(0);
          }
        }

        // Nodes
        ibuf[0] = FV_NODES;
        ibuf[1] = nn;
        fwrite(ibuf,isize,2,fp);
        ftmp = (float*)malloc(ibuf[1]*fsize);
        // X coordinate
        for (n=0; n < ibuf[1]; n++)
          ftmp[n] = (float)(node[n].vert[0]);
        fwrite(ftmp,fsize,ibuf[1],fp);
        // Y coordinate
        for (n=0; n < ibuf[1]; n++)
          ftmp[n] = (float)(node[n].vert[1]);
        fwrite(ftmp,fsize,ibuf[1],fp);
        // Z coordinate
        for (n=0; n < ibuf[1]; n++)
          ftmp[n] = (float)(node[n].vert[2]);
        fwrite(ftmp,fsize,ibuf[1],fp);

        free(ftmp);

        // Boundary faces
        for (b=0; b < nb; b++)
        {
          nt = nq = 0;
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = abs(boundary[b].polygon_list->list[i]);
            if (face[f].node_list->max == 3)
              nt++;
            else if (face[f].node_list->max == 4)
              nq++;
          }
          ibuf[0] = FV_FACES;
          ibuf[1] = b+1;
          ibuf[2] = nt+nq;
          fwrite(ibuf,isize,3,fp);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            k = face[abs(f)].node_list->max;
            if (k < 3 || k > 4) continue;
            n0 = face[abs(f)].node_list->list[0]+1;
            n1 = face[abs(f)].node_list->list[1]+1;
            n2 = face[abs(f)].node_list->list[2]+1;
            if (k == 4)
              n3 = face[abs(f)].node_list->list[3]+1;
            else
              n3 = 0;
            if (f < 0)
            {
              ibuf[0] = n2;
              ibuf[1] = n1;
              ibuf[2] = n0;
              ibuf[3] = n3;
            } else
            {
              ibuf[0] = n0;
              ibuf[1] = n1;
              ibuf[2] = n2;
              ibuf[3] = n3;
            }
            fwrite(ibuf,isize,4,fp);
          }
        }
        for (b=0; b < nb; b++)
        {
          nt = nq = 0;
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = abs(boundary[b].polygon_list->list[i]);
            if (face[f].node_list->max == 3)
              nt++;
            else if (face[f].node_list->max == 4)
              nq++;
          }
          if (nt+nq == boundary[b].polygon_list->max)
            continue;
          ibuf[0] = FV_ARB_POLY_FACES;
          ibuf[1] = b+1;
          ibuf[2] = boundary[b].polygon_list->max-nt-nq;
          fwrite(ibuf,isize,3,fp);
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            k = face[abs(f)].node_list->max;
            if (k == 3 || k == 4) continue;
            k=0;
            ibuf[k++] = face[abs(f)].node_list->max;
            if (f < 0)
              for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                ibuf[k++] = face[abs(f)].node_list->list[j]+1;
            else
              for (j=0; j < face[f].node_list->max; j++)
                ibuf[k++] = face[f].node_list->list[j]+1;
            fwrite(ibuf,isize,k,fp);
          }
        }

        // Element section
        int walls[6];
        walls[0] = NOT_A_WALL;
        walls[1] = NOT_A_WALL;
        walls[2] = NOT_A_WALL;
        walls[3] = NOT_A_WALL;
        walls[4] = NOT_A_WALL;
        walls[5] = NOT_A_WALL;
        ntet = npyr = npri = nhex = nply = 0;
        for (e=0; e < nelems; e++)
        {
          switch(element_type(e))
          {
            case 4: ntet++; break;
            case 5: npyr++; break;
            case 6: npri++; break;
            case 8: nhex++; break;
            default: nply++; break;
          }
        }
        ibuf[0] = FV_ELEMENTS;
        ibuf[1] = ntet; // # of tetrahedron
        ibuf[2] = nhex; // # of hexahedron
        ibuf[3] = npri; // # of prism
        ibuf[4] = npyr; // # of pyramid
        fwrite(ibuf,isize,5,fp);
        for (e=0; e < nelems; e++)
        {
          switch(element_type(e))
          {
            case 4:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              ibuf[0]=fv_encode_elem_header(FV_TET_ELEM_ID,walls);
              fwrite(ibuf,isize,1,fp);
              ibuf[0]=n0;
              ibuf[1]=n2;
              ibuf[2]=n1;
              ibuf[3]=n3;
              fwrite(ibuf,isize,4,fp);
              break;
            case 5:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              ibuf[0]=fv_encode_elem_header(FV_PYRA_ELEM_ID,walls);
              fwrite(ibuf,isize,1,fp);
              ibuf[0]=n0;
              ibuf[1]=n1;
              ibuf[2]=n2;
              ibuf[3]=n3;
              ibuf[4]=n4;
              fwrite(ibuf,isize,5,fp);
              break;
            case 6:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              n5 = element[e].node_list->list[5]+1;
              ibuf[0]=fv_encode_elem_header(FV_PRISM_ELEM_ID,walls);
              fwrite(ibuf,isize,1,fp);
              ibuf[0]=n0;
              ibuf[1]=n3;
              ibuf[2]=n5;
              ibuf[3]=n1;
              ibuf[4]=n2;
              ibuf[5]=n4;
              fwrite(ibuf,isize,6,fp);
              break;
            case 8:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              n5 = element[e].node_list->list[5]+1;
              n6 = element[e].node_list->list[6]+1;
              n7 = element[e].node_list->list[7]+1;
              ibuf[0]=fv_encode_elem_header(FV_HEX_ELEM_ID,walls);
              fwrite(ibuf,isize,1,fp);
              ibuf[0]=n0;
              ibuf[1]=n1;
              ibuf[2]=n3;
              ibuf[3]=n2;
              ibuf[4]=n4;
              ibuf[5]=n5;
              ibuf[6]=n7;
              ibuf[7]=n6;
              fwrite(ibuf,isize,8,fp);
              break;
            default:
              break;
          }
        }

        if (nply > 0)
        {
          ibuf[0] = FV_ARB_POLY_ELEMENTS;
          ibuf[1] = nply;
          fwrite(ibuf,isize,2,fp);
          for (e=0; e < nelems; e++)
          {
            k= element_type(e);
            if (k==4 || k == 5 || k == 6 || k == 8) continue;
            ibuf[0] = element[e].face_list->max;  // number of faces
            ibuf[1] = element[e].node_list->max; // number of nodes in element
            ibuf[2] = -1;  // center node
            fwrite(ibuf,isize,3,fp);
            for (i=0; i < element[e].face_list->max; i++)
            {
              f = element[e].face_list->list[i];
              ibuf[0] = NOT_A_WALL;
              fwrite(ibuf,isize,1,fp);  // wall value
              k=0;
              ibuf[k++] = face[abs(f)].node_list->max;  // number of face nodes
              if (f < 0)
                for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                  ibuf[k++] = face[abs(f)].node_list->list[j]+1;  // face nodes
              else
                for (j=0; j < face[f].node_list->max; j++)
                  ibuf[k++] = face[f].node_list->list[j]+1;  // face nodes
              fwrite(ibuf,isize,k,fp);
              ibuf[0] = 0;
              fwrite(ibuf,isize,1,fp); // number of hanging nodes
            }
          }
        }

        ftmp = (float*)malloc(nn*fsize);
        ibuf[0] = FV_VARIABLES;
        fwrite(ibuf,isize,1,fp);
        for (i=0; i < num_func; i++)
        {
          for (n=0; n < nn; n++)
            ftmp[n] = func[i][n];
          fwrite(ftmp,fsize,nn,fp);
        }

        free(ftmp);

        ibuf[0] = FV_BNDRY_VARS;
        fwrite(ibuf,isize,1,fp);

        fclose(fp);
      }

      break;
    case 2: // Write ASCII
      sprintf(filename,"%s.crunch",buff);
      for (proc = 0; proc < num_procs; proc++)
      {
        #ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        if (proc != my_rank)
          continue;

        // Open file for write
        if (proc == 0)
        {
          if ((fp = fopen(filename,"w")) == 0)
          {
            fprintf(out_f,"\nError opening file <%s>.",filename);
            #ifdef PARALLEL
            MPI_Abort(MPI_COMM_WORLD,0);
            #endif
            exit(0);
          }
          // name and version number
          fprintf(fp,"FieldView_Grids 3 0\n");

          fprintf(fp,"Constants\n");
          ftmp = (float*)malloc(4*fsize);
          ftmp[0] = 0.0; // Time
          ftmp[1] = 0.0; // Mach
          ftmp[2] = 0.0; // Alpha
          ftmp[3] = 0.0; // Reynold's number
          fprintf(fp,"%g %g %g %g\n",ftmp[0],ftmp[1],ftmp[2],ftmp[3]);
          free(ftmp);

          // Number of grids
          fprintf(fp,"Grids\n");
          fprintf(fp,"%d\n",num_procs);

          // Number of boundary types
          fprintf(fp,"Boundary Table\n%d\n",nb);

          // Boundary types
          for (b=0; b < nb; b++)
          {
            for (i=0; i < 80; i++)
              buff[i] = ' ';
            sprintf(buff,"%s",boundary[b].name);
            fprintf(fp,"1 0 1 %s\n",buff);
          }

          // Number of variables
          fprintf(fp,"Variable Names\n%d\n",num_func);
          for (i=0; i < num_func; i++)
            fprintf(fp,"%s\n",func_name[i]);
          fprintf(fp,"Boundary Variable Names\n0\n");

        } else
        {
          if ((fp = fopen(filename,"a")) == 0)
          {
            fprintf(out_f,"\nError opening file <%s>.",filename);
            #ifdef PARALLEL
            MPI_Abort(MPI_COMM_WORLD,0);
            #endif
            exit(0);
          }
        }

        // Node information  
        fprintf(fp,"Nodes\n");
        fprintf(fp,"%d\n",nn);
        for (n=0; n < nn; n++)
          fprintf(fp,"%22.15e %22.15e %22.15e\n",node[n].vert[0],node[n].vert[1],node[n].vert[2]);

        // Boundary section
        n=0;
        for (b=0; b < nb; b++)
          n += boundary[b].polygon_list->max;
        fprintf(fp,"Boundary Faces\n");
        fprintf(fp,"%d\n",n);
        for (b=0; b < nb; b++)
        {
          for (i=0; i < boundary[b].polygon_list->max; i++)
          {
            f = boundary[b].polygon_list->list[i];
            fprintf(fp,"%d %d",b+1,face[abs(f)].node_list->max);
            if (f < 0)
            {
              for (j = face[abs(f)].node_list->max-1; j >= 0; j--)
                fprintf(fp," %d",face[abs(f)].node_list->list[j]+1);
            } else
            {
              for (j=0; j < face[f].node_list->max; j++)
                fprintf(fp," %d",face[f].node_list->list[j]+1);
            }
            fprintf(fp,"\n");
          }
        }

        // Element section
        fprintf(fp,"Elements\n");
        for (e=0; e < nelems; e++)
        {
          switch(element_type(e))
          {
            case 4:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              fprintf(fp,"1 1\n%d %d %d %d\n",n0,n2,n1,n3);
              break;
            case 5:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              fprintf(fp,"4 1\n%d %d %d %d %d\n",n0,n1,n2,n3,n4);
              break;
            case 6:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              n5 = element[e].node_list->list[5]+1;
              fprintf(fp,"3 1\n%d %d %d %d %d %d\n",n0,n3,n5,n1,n2,n4);
              break;
            case 8:
              n0 = element[e].node_list->list[0]+1;
              n1 = element[e].node_list->list[1]+1;
              n2 = element[e].node_list->list[2]+1;
              n3 = element[e].node_list->list[3]+1;
              n4 = element[e].node_list->list[4]+1;
              n5 = element[e].node_list->list[5]+1;
              n6 = element[e].node_list->list[6]+1;
              n7 = element[e].node_list->list[7]+1;
              fprintf(fp,"2 1\n%d %d %d %d %d %d %d %d\n",n0,n1,n3,n2,n4,n5,n7,n6);
              break;
            default:
              fprintf(fp,"5 1\n");
              fprintf(fp,"%d %d -1\n",element[e].face_list->max,element[e].node_list->max);
              for (i=0; i < element[e].face_list->max; i++)
              {
                f = element[e].face_list->list[i];
                fprintf(fp,"%d ",face[abs(f)].node_list->max);
                if (f < 0)
                {
                  for (j=face[abs(f)].node_list->max-1; j >= 0; j--)
                    fprintf(fp,"%d ",face[abs(f)].node_list->list[j]+1);
                } else
                {
                  for (j=0; j < face[f].node_list->max; j++)
                    fprintf(fp,"%d ",face[f].node_list->list[j]+1);
                }
                fprintf(fp,"0\n");
              }
              break;
          }
        }

        // Variable section
        fprintf(fp,"Variables\n");

        for (i=0; i < num_func; i++)
          for (n=0; n < nn; n++)
            fprintf(fp,"%22.15e\n",func[i][n]);

        // Boundary variable section
        fprintf(fp,"Boundary Variables\n");

        fclose(fp);
      }

      break;
    default:
      break;
  }

  return;

}
