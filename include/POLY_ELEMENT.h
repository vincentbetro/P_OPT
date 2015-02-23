#include <stdio.h>
#include "Pmap.h"
#include "List.h"

#ifndef POLY_ELEMENT_h
#define POLY_ELEMENT_h

class POLYGONAL_ELEMENT
{
 public:
  POLYGONAL_ELEMENT() // Constructor
  {
    node_list = 0;
    return;
  }
  ~POLYGONAL_ELEMENT() // Destructor
  {
    if (node_list != 0)
      delete node_list;
    node_list = 0;
    return;
  }
  void Add_Nodes(List *nlist)
  {
    if (node_list == 0)
      node_list = new List();
    *node_list = *nlist;
  }
  void Initialize()
  {
    node_list = 0;
    return;
  }
  void Destroy()
  {
    if (node_list != 0)
      delete node_list;
    node_list = 0;
    return;
  }

  List *node_list; // list of nodes in face
};

class POLYHEDRAL_ELEMENT
{
 public:
  POLYHEDRAL_ELEMENT() // Constructor
  {
    face_list = new List();
    node_list = new List();
    me = Pmap(-1,-1,-1);
    return;
  }
  ~POLYHEDRAL_ELEMENT() // Destructor
  {
    if (face_list != 0)
      delete face_list;
    face_list = 0;
    if (node_list != 0)
      delete node_list;
    node_list = 0;
    me = Pmap(-1,-1,-1);
    return;
  }
  void Destroy()
  {
    if (face_list != 0)
      delete face_list;
    face_list = 0;
    if (node_list != 0)
      delete node_list;
    node_list = 0;
    me = Pmap(-1,-1,-1);
    return;
  }
  void Initialize() // for Resetting or intializing an array of elements.
  {
    face_list = new List();
    node_list = new List();
    me = Pmap(-1,-1,-1);
    return;
  }
  void Add_Face(int fnum)
  {
    if (face_list == 0)
      face_list = new List();
    face_list->Add_To_List(fnum);
    return;
  }

  List *face_list; // list of faces
  List *node_list; // list of nodes
  Pmap me;    // processor map
};

#endif
