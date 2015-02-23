#include <stdio.h>

#ifndef Point_h
#define Point_h
class Point
{
  public:
    inline Point(double x=0.0, double y=0.0, double z=0.0)
    { pos[0] = x; pos[1] = y; pos[2] = z; }
  ~Point() { }
  void print(FILE *outf)
  {
    fprintf(outf,"\nPoint (x,y,z)= (%.12g,%.12g,%.12g)",pos[0],pos[1],pos[2]);
  }
  inline Point operator + (const Point &) const;
  inline Point operator - (const Point &) const;
  inline Point &operator += (const Point &);
  inline Point &operator -= (const Point &);
  inline Point &operator += (double);
  inline Point &operator -= (double);
  inline Point &operator *= (double);
  inline Point &operator /= (double);
  inline Point operator * (double) const;
  inline Point operator / (double) const;
  inline double & operator () (int);
  inline double operator () (int) const;
  inline double & operator [] (int);

  double pos[3];
};

inline double & Point::operator () (int i) 
{ return pos[i]; }

inline double Point::operator () (int i) const
{ return pos[i]; }

inline double & Point::operator [] (int i) 
{ return pos[i]; }

inline Point Point::operator * (double other) const 
{ return Point(pos[0]*other, pos[1]*other, pos[2]*other); }

inline Point Point::operator / (double other) const 
{
#ifdef _DEBUG_
	if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
  return Point(pos[0]/other, pos[1]/other, pos[2]/other);
} 

inline Point Point::operator + (const Point &other) const 
{ return Point(pos[0]+other.pos[0], pos[1]+other.pos[1], pos[2]+other.pos[2]); }

inline Point Point::operator - (const Point &other) const 
{ return Point(pos[0]-other.pos[0], pos[1]-other.pos[1], pos[2]-other.pos[2]); }

inline Point& Point::operator += (const Point &other) 
{
  for (int i = 0; i < 3; i++)
    pos[i]+=other.pos[i];
  return *this;
}

inline Point& Point::operator -= (const Point &other) 
{
  for (int i = 0; i < 3; i++)
    pos[i]-=other.pos[i];
  return *this;
}

inline Point& Point::operator += (double other) 
{
  for (int i = 0; i < 3; i++)
    pos[i]+=other;
  return *this;
}

inline Point& Point::operator -= (double other) 
{
  for (int i = 0; i < 3; i++)
    pos[i]-=other;
  return *this;
}

inline Point& Point::operator *= (double other) 
{
  for (int i = 0; i < 3; i++)
    pos[i]*=other;
  return *this;
}

inline Point& Point::operator /= (double other) 
{
#ifdef _DEBUG_
	if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
  for (int i = 0; i < 3; i++)
    pos[i]/=other;
  return *this;
}

#endif
