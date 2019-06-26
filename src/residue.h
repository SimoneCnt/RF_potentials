#ifndef _RESIDUE_H_
#define _RESIDUE_H_

#include <string>
#include <map>
using namespace std;

#include "atom.h"

#include "chain.h"

class Residue {
public:
  Chain *chain;		// our chain
  Residue *prev;	//previous residue
  Residue *next;	//next residue
  char chain_name;	//chain symbol
  int model;		// model number
  char res;		// residue name in 1-letter code
  char *name;		// residue name in 3-letter code
  int res_idx;		// numerical residue index (for compartability with Chain)
  int number;		//residue_number from pdb file
  char icode;		//insertion iCODE

  int n_atoms;		//number of atoms of residue
  Atom *atoms[NATM_MAX]; // list of pointers to the residue Atom atoms
  map< string, Atom *> atom; // named list of pointers to the residue Atom atoms

  Residue( Chain *C, int resno );
  ~Residue(void);
  void MakeVirtualCb( double K = 1 );
  void WritePDB( int i, int &i_seq, int &mdl, FILE *out );
};
#endif


