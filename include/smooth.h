#include <stdio.h>
#include "Pmap.h"
#include "Point.h"
#include "Vector.h"
#include "geometry.h"
#include "List.h"
#include "Bnormal.h"
#include "POLY_ELEMENT.h"
#include "NODE.h"

#ifndef Poly_mesh_objects_h
#define Poly_mesh_objects_h

class SNODE: public NODE
{
 public:
  SNODE()
  {
    me = Pmap(-1,-1,-1);
    nabor = 0;
    element = 0;
    return;
  }
  ~SNODE()
  {
    me = Pmap(-1,-1,-1);
    if (nabor != 0) delete nabor;
    if (element != 0) delete element;
    return;
  }
  void print(FILE *prnt)
  {
    vert.print(prnt);
    me.print(prnt);
  }

  List *nabor;     // List of neighboring nodes
  List *element;   // List of surrounding elements
};

class BOUNDARY_MESH
{
 public:
  BOUNDARY_MESH()
  {
    polygon_list=0;
    name=NULL;
    return;
  }
  ~BOUNDARY_MESH()
  {
    if (polygon_list)
      delete polygon_list;
    if (name != NULL)
      free(name);
    polygon_list=0;
    name=NULL;
    return;
  }
  void Initialize() // for Resetting or intializing an array of boundaries.
  {
    polygon_list=0;
    name=NULL;
    return;
  }
  void Destroy()
  {
    if (polygon_list)
      delete polygon_list;
    if (name != NULL)
      free(name);
    polygon_list=0;
    name=NULL;
    return;
  }

  List *polygon_list;
  char *name;
};

class POLYMESH
{
 public:
  POLYMESH()
  {
    nn=nelems=nfaces=nb=0;
    node_dim=element_dim=face_dim=boundary_dim=0;
    node=0;
    element=0;
    face=0;
    boundary=0;
    return;
  }
  ~POLYMESH()
  {
    if (node > 0) free(node);
    if (element > 0)
    {
      for (int e=0; e < nelems; e++)
        element[e].Destroy();
      free(element);
    }
    if (face > 0)
    {
      for (int f=0; f <= nfaces; f++)
        face[f].Destroy();
      free(face);
    }
    if (boundary > 0)
    {
      for (int b=0; b < nb; b++)
        boundary[b].Destroy();
      free(boundary);
    }
    nn=nelems=nfaces=nb=0;
    node_dim=element_dim=face_dim=boundary_dim=0;
    node=0;
    element=0;
    face=0;
    boundary=0;
    return;
  }
  int check_volumes();
  int check_Jacobians();
  void create_element_maps();
  int critical_points(geometry *geom, int **cpn);
  double element_volume(int e, Point &cg);
  void mesh_stats();
  bool node_in_element(int n, int e);
  int read_mesh(char sname[]);
  int write_mesh(char sname[]);
  int smooth_io(int mode, geometry *geom, char sname[]);
  void exchange_node_double(double *array);
  void exchange_node_int(int *array);
  void exchange_node_vector(Vector *array);
  int adjacent_nodes(int e, int n, List *nlist, List *flist);
  int element_type(int e);
  int face_number(List *nlist, List **flist);
  void Fieldview(int mode, char sname[], int num_func, char (*func_name)[80], double **func);
  void Ensight(int mode, char sname[], int num_func, char (*func_name)[80], double **func,
               int nl=0, int *ne=0, int ***edge=0);
  double cost_function(int n);
  double cost_function(geometry *geom, int n, List *flist, int bindex[], List **blist);
  int opt_smooth(geometry *geom, int nsmoo, int mvbnd, int mvint, int sflag, int eflag, double threshold);

  int nn, nelems, nfaces, nb;
  int node_dim, element_dim, face_dim, boundary_dim;
  SNODE *node;
  POLYHEDRAL_ELEMENT *element;
  POLYGONAL_ELEMENT *face;
  BOUNDARY_MESH *boundary;
};
#endif
