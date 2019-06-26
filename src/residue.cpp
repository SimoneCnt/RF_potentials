#include <stdlib.h>
#include "residue.h"
#include "point_3.h"


Residue::Residue( Chain *pC, int resno ) {
  int i, j;

  for( j = 0; j < NATM_MAX; j++ )
    atoms[j] = NULL;
  prev = next = NULL;
  atom["CB"] = atom["CD"] = atom["CD1"] = atom["CD2"] = atom["CE"] = atom["CE1"] = NULL;
  atom["CE2"] = atom["CE3"] = atom["CG"] = atom["CG1"] = atom["CG2"] = atom["CH2"] = NULL;
  atom["CZ"] = atom["CZ2"] = atom["CZ3"] = atom["ND1"] = atom["ND2"] = atom["NE"] = NULL;
  atom["NE1"] = atom["NE2"] = atom["NH1"] = atom["NH2"] = atom["NZ"] = atom["OD1"] = NULL;
  atom["OD2"] = atom["OE1"] = atom["OE2"] = atom["OG"] = atom["OG1"] = atom["OH"] = NULL;
  atom["SD"] = atom["SG"] = NULL;

  chain = pC;
  chain_name = pC->residue_chain[resno];
  res = pC->sequence[resno];
  name = pC->residue_name[resno];
  number = pC->residue_number[resno];
  icode = pC->residue_icode[resno];
  model = pC->residue_model[resno];

  res_idx = pC->iresidues[resno];
  n_atoms = pC->aa[res_idx].size;

  for( j = 0; j < n_atoms; j++ ) {
    i = pC->ica[resno] + j;
    atoms[j] = new Atom( pC, this, pC->aa[res_idx].atoms[j], i, 
			 pC->pts[0][i], pC->pts[1][i], pC->pts[2][i], pC->tfactor[i], pC->faked_atom[i] );
  }

  atom["N"] = atoms[0];
  atom["CA"] = atoms[1];
  atom["C"] = atoms[2];
  atom["O"] = atoms[3];
  if( res != 'G' )
    atom["CB"] = atoms[4];
  switch( res ) {
  case 'C' :
    atom["SG"] = atoms[5];
    break;
  case 'S' :
    atom["OG"] = atoms[5];
    break;
  case 'T' :
    atom["OG1"] = atoms[5];
    atom["CG2"] = atoms[6];
    break;
  case 'P' :
    atom["CG"] = atoms[5];
    atom["CD"] = atoms[6];
    break;
  case 'V' :
    atom["CG1"] = atoms[5];
    atom["CG2"] = atoms[6];
    break;
  case 'L' :
    atom["CG"] = atoms[5];
    atom["CD1"] = atoms[6];
    atom["CD2"] = atoms[7];
    break;
  case 'I' :
    atom["CG1"] = atoms[5];
    atom["CG2"] = atoms[6];
    atom["CD1"] = atoms[7];
    break;
  case 'N' :
    atom["CG"] = atoms[5];
    atom["OD1"] = atoms[6];
    atom["ND2"] = atoms[7];
    break;
  case 'D' :
    atom["CG"] = atoms[5];
    atom["OD1"] = atoms[6];
    atom["OD2"] = atoms[7];
    break;
  case 'M' :
    atom["CG"] = atoms[5];
    atom["SD"] = atoms[6];
    atom["CE"] = atoms[7];
    break;
  case 'E' :
    atom["CG"] = atoms[5];
    atom["CD"] = atoms[6];
    atom["OE1"] = atoms[7];
    atom["OE2"] = atoms[8];
    break;
  case 'Q' :
    atom["CG"] = atoms[5];
    atom["CD"] = atoms[6];
    atom["OE1"] = atoms[7];
    atom["NE2"] = atoms[8];
    break;
  case 'K' :
    atom["CG"] = atoms[5];
    atom["CD"] = atoms[6];
    atom["CE"] = atoms[7];
    atom["NZ"] = atoms[8];
    break;
  case 'H' :
    atom["CG"] = atoms[5];
    atom["ND1"] = atoms[6];
    atom["CD2"] = atoms[7];
    atom["CE1"] = atoms[8];
    atom["NE2"] = atoms[9];
    break;
  case 'F' :
    atom["CG"] = atoms[5];
    atom["CD1"] = atoms[6];
    atom["CD2"] = atoms[7];
    atom["CE1"] = atoms[8];
    atom["CE2"] = atoms[9];
    atom["CZ"] = atoms[10];
    break;
  case 'R' :
    atom["CG"] = atoms[5];
    atom["CD"] = atoms[6];
    atom["NE"] = atoms[7];
    atom["CZ"] = atoms[8];
    atom["NH1"] = atoms[9];
    atom["NH2"] = atoms[10];
    break;
  case 'Y' :
    atom["CG"] = atoms[5];
    atom["CD1"] = atoms[6];
    atom["CD2"] = atoms[7];
    atom["CE1"] = atoms[8];
    atom["CE2"] = atoms[9];
    atom["CZ"] = atoms[10];
    atom["OH"] = atoms[11];
    break;
  case 'W' :
    atom["CG"] = atoms[5];
    atom["CD1"] = atoms[6];
    atom["CD2"] = atoms[7];
    atom["NE1"] = atoms[8];
    atom["CE2"] = atoms[9];
    atom["CE3"] = atoms[10];
    atom["CZ2"] = atoms[11];
    atom["CZ3"] = atoms[12];
    atom["CH2"] = atoms[13];
    break;
  }
}

void Residue::MakeVirtualCb( double K ) {
  if( res_idx != 7 && atoms[4]->faked == false )             // this function is for GLY only!
    return;
  Point_3 Ca( atoms[1]->x[0], atoms[1]->x[1], atoms[1]->x[2] );
  Point_3 p0( (atoms[1]->x[0] - atoms[0]->x[0]), (atoms[1]->x[1] - atoms[0]->x[1]), 
	      (atoms[1]->x[2] - atoms[0]->x[2]) );
  Point_3 p1( (atoms[2]->x[0] - atoms[1]->x[0]), (atoms[2]->x[1] - atoms[1]->x[1]), 
	      (atoms[2]->x[2] - atoms[1]->x[2]) );
  Point_3 p01 = p1 - p0;
  Point_3 p1n = (p0 >> p1 );

  p01 = p01 / ~p01;
  p1n = p1n / ~p1n;
  Point_3 p2n = (p01 >> p1n );
  p2n = p2n / ~p2n;

//    Point_3 Cb = - p01*0.8833459118601273 - p1n*1.2492397688194208; // ideal. does not work.
  Point_3 Cb = p2n*0.029 - p01*0.940 - p1n*1.207; // optimized: mean deviation is 0.078924 +/- 0.058336;
  Cb = Cb * K;
  Cb = Cb + Ca;
  if( res_idx == 7 ) {
    atoms[4] = new Atom( chain, this, (char*)"CB", -1, Cb.x, Cb.y, Cb.z, 999., true );
    atom["CB"] = atoms[4];
  }
  else {
    atoms[4]->x[0] = Cb.x;
    atoms[4]->x[1] = Cb.y;
    atoms[4]->x[2] = Cb.z;
    atoms[4]->faked = false;
  }
}
		
Residue::~Residue(void) {
  for( int j = 0; j < NATM_MAX; j++ )
    delete atoms[j];
  atom.clear();
}

void Residue::WritePDB( int i, int &i_seq, int &mdl, FILE *out ) {
  int j;
  char tmp[255];
  int n_models = chain->n_models;
  int n_residues = chain->n_residues;

  if( n_models > 1 )
    if( mdl != model ) {
      mdl = model;
      fprintf( out, "MODEL      %3d\n", mdl+1 );
    }

  for( j = 0; j < n_atoms; j++ ) {
    i_seq++;
    sprintf(tmp, "ATOM %6d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f  1.00%6.2f",
	    i_seq, atoms[j]->name, name, chain_name, number,
	    atoms[j]->x[0], atoms[j]->x[1], atoms[j]->x[2], atoms[j]->tfactor );
    fprintf(out, "%s\n", tmp);
  }
  if( res_idx == 7 && atoms[4] != NULL ) {
    i_seq++;
    sprintf(tmp, "ATOM %6d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f  1.00%6.2f",
	    i_seq, atoms[4]->name, name, chain_name, number,
	    atoms[4]->x[0], atoms[4]->x[1], atoms[4]->x[2], atoms[4]->tfactor );
    fprintf(out, "%s\n", tmp);
  }

  if( n_models > 1 )
    if( i+1 < n_residues )
      if( mdl != next->model )
	fprintf( out, "ENDMDL\n" );
  i++;

//0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//ATOM    157  CD2 LEU A  19       1.640  27.834  22.650  1.00 13.64           C  
  if( n_models > 1 )
    fprintf( out, "ENDMDL\n" );
}
