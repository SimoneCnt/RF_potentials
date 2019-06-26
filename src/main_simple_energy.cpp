#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "residue.h"
#include "energy.h"


int main(int argc, char *argv[]) {
  int i;
  FILE *f;
  char buff[1024];

  if( argc != 3 ) {
    printf( "2 args expected, %d found:\n", argc-1 );
    for( i = 1; i < argc; i++ )
      printf( "\t%d\t\"%s\"\n", i, argv[i] );
    printf( "\nUsage: %s <PDB_file> <potential_file>\n", argv[0] );
    exit( 1 );
  }

// create Potential object array
  int np, nbin, Npoten = 1;

  Chain *pdb = new Chain( '*', argv[1] );
  i = pdb->ReadPDB();
  if( i == 0 ) {
    printf( "Can't read PDB file %s\n", argv[1] );
    exit( 1 );
  }
   
  Energy *e = new Energy( argv[2] );
  Residue **r1 = pdb->MkRes();
  
  printf( "%g\n", e->Execute( r1 ) );
  delete e;
  for( i = 0; i < pdb->n_residues; i++ )
    delete r1[i];
  delete r1;
  delete pdb;

  return 0;
}
