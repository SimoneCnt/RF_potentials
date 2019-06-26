#ifndef _CHAIN_
#define _CHAIN_

#include <vector>
using namespace std;

#include "aminoacid.h"
#include "residue.h"

class Chain : public Aminoacid {

  vector<double> tmp_tfactor;
  vector<double> tmp_pts[3];
  vector<int> a_idx;
  bool find_atm[NATM_MAX];
  int FinalizeResidue( int &i_atom, int &i_residue, int &loaded_atoms );
 public:
  int n_aa;			// number of aminoacids
  int n_models;			// number of models in the protein
  int n_residues;		// number of complete residues of protein
  double resolution;		// X-ray resolution. == -1 for NMR entries

  int short_residues;		// number of incomplete residues
  int n_atoms;			// number of atoms of protein
  char *pdb_file;		// protein file
  char pdb_code[5];		// pdb code
  char chain_symbol;		// chain name
  vector<char> sequence;	// aminoacid sequence [0 - n_residues]
  vector<char> residue_chain;	// chain for each residue
  vector<int> residue_model;	// model number for each residue
  vector<int> residue_number;	// residue_number from pdb file [0 - n_residues]
  vector<int> short_list;	// residue_number from pdb file [0 - n_residues]
  vector<char*> residue_name;	// residue name in 3-letter code
  vector<char> residue_icode;	// insertion iCODE
  vector<int> ica;		// index of first atom of residue in pts [0 - n_residues]
  vector<int> iresidues;	// index of residue [0-19] 0 - ALA, 1 - ARG, ... (see aminoacid.cpp) [0 - n_residues]
  vector<int> faked_atom;	// list of faked atoms in the incomplete residues [0 - n_atoms]
  vector<double> pts[3];	// coordinates of all atoms of protein x=pts[0][i], y=pts[1][i], z=pts[2][i]; i = [0 - n_atoms]
  vector<double> tfactor;	// temperature factor [0 - n_atoms]

  static Aminoacid aa[21];		//amino acid array

  Chain( Chain *c );
  Chain( char symbol, char *file );
  ~Chain(void);
  int ReadPDB(void);
  int WritePDB(char *filename);

  Residue **MkRes(void);
};
#endif
