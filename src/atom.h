#ifndef _ATOM_H_
#define _ATOM_H_

#include <math.h>

class Chain;
class Residue;

class Atom {
public:
  Chain *chain; // our chain
  Residue *res; // our residue

  char *name;	// atom name
  int no;	// atom_number from pdb file
  double x[3];  // coordinates
  double tfactor;//temperature factor
  bool faked;

  Atom(Chain *C, Residue *R, char *Name, int No, double X, double Y, double Z, double T, bool F ) { 
    chain = C;
    res = R;
    name = Name;
    no = No;
    x[0] = X; x[1] = Y; x[2] = Z;  
    tfactor = T;
    faked = F;
  }
  ~Atom(void) {};

  double Dist2(Atom *a) {
    double xx = a->x[0] - x[0];
    double yy = a->x[1] - x[1];
    double zz = a->x[2] - x[2];
    return xx*xx+yy*yy+zz*zz;
  }

  double Distance(Atom *a) {
    return sqrt( Dist2(a) );
  }

};
#endif
