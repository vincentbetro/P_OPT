#include <stdio.h>
#include "Pmap.h"
#include "Point.h"

#ifndef NODE_h
#define NODE_h

class NODE
{
 public:
  NODE()
  {
    me = Pmap(-1,-1,-1);
    return;
  }
  ~NODE()
  {
    me = Pmap(-1,-1,-1);
    return;
  }
  void print(FILE *prnt)
  {
    vert.print(prnt);
    me.print(prnt);
  }

  Point vert;      // Physical coordinate triplet
  Pmap me;         // Parallel map info
};

#endif
