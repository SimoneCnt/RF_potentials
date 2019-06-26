#ifndef _POINT3_H_
#define _POINT3_H_

#include <math.h>

class Point_3 { 
 public:
  double x,y,z;
  Point_3(void) {x = y = z = 0.;}
  Point_3(double dx, double dy, double dz) {x=dx; y=dy; z=dz;}

  Point_3 operator/(double norm) {x=x/norm; y=y/norm; z=z/norm; return *this;}
  Point_3 operator*(double norm) {x=x*norm; y=y*norm; z=z*norm; return *this;}
  Point_3 operator+(Point_3 p) {return Point_3(x+p.x, y+p.y, z+p.z);}
  Point_3 operator-(Point_3 p) {return Point_3(x-p.x, y-p.y, z-p.z);}
  double operator*(Point_3 p) {return (x*p.x + y*p.y + z*p.z);}

  Point_3 operator/=(double norm) {x=x/norm; y=y/norm; z=z/norm; return *this;}
  double operator~() {return sqrt(x*x+y*y+z*z);}
  Point_3 operator>>(Point_3 p) {return Point_3(y*p.z - z*p.y, -x*p.z + z*p.x, x*p.y - y*p.x);}

};
#endif
