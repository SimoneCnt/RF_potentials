#ifndef _ENERGY_H_
#define _ENERGY_H_
#ifndef PI
#define PI 3.1415926
#endif

#define MAXATMTYPES 5
#define NATM 167  // maximal number of possible atom types

#define RT 1
#define MY_INF 10.

#include <openssl/evp.h>
#include <math.h>

class Residue;

class Energy {
//  POTENTIAL properties
  int na[MAXATMTYPES];
  int Nprot;            // number of proteins in the statistics
  bool Processed;       // Final calculation completeness indicator
  double BinSize;

  unsigned char md_value[EVP_MAX_MD_SIZE]; // openSSL sha1 digest for the poten to ensure I/O integrity

  double aa_ttl;           // total number of the amino acid residues in the potential statistics
  double pair_total;       // total number of the all residue pairs over all distances used to derive potential

  int AverType;

// PROTEIN properties
  Chain *pdb;
  int n_atoms;
  int Nres;
  Residue **r;

  bool checksum_verify(void);
  unsigned int checksum(unsigned char *md_value);
  void PrintSignature( FILE *out = stdout );
  void PrintSignature( FILE *out, unsigned char *md_val );

  inline int BinMin( double DistMin ) { return ( DistMin == -1. ) ? 0 : (int)(( DistMin - MinDist) / (MaxDist - MinDist) * Nbins); }
  inline int BinMax( double DistMax ) { return ( DistMax == -1. ) ? 0 : (int)(( DistMax - MinDist) / (MaxDist - MinDist) * Nbins); }
  double SumUp( double DistMin, double DistMax );


 public:
//  POTENTIAL properties
  char Descr[128]; // brief potential description
  int NatmType;	   // type of interaction centers: 0="CA", 1="AllHeavyAtoms", 2="unused", 3="unused", 4="CB", 5="CB_noG", 6="unused", 7="unused",

  bool Directed;   // orientation dependence property
  int natm;	   // number of interaction centes ("atoms")

  int ChainSepMin; // minimal chain separation
  int ChainSepMax; // maximal chain separation

  double MinDist;  // minimal distance cutoff for the pair to be taken into account
  double AverDist; // unused
  double MaxDist;  // maximal distance cutoff
  int Nbins;       // number of space bins
  double AppCutoff;// maximal distance of the potential applicability
  int ChainDependent;   // Potential type: {0,1} - whole protein; {2,3} - chain only; {0,2} - symmentic; {1,3} - assymetric;
  double ***poten; // potential matrix
  int DepthDep;    // unused
  bool Smooth;     // unused

// PROTEIN energy
  double **energy;  // energy calculated for each interaction center and bin
  double E;	    // total energy

  Energy( char *PotenFile );
  Energy( const Energy &e );

  double Execute( Residue **R );

  bool ReadPotential( const char *FileName );
  void PrintPotential( FILE *out, int mode = 0 );

  ~Energy(void);
};
#endif
