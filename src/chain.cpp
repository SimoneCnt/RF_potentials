#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef __unix__
#include <zlib.h>
#else 
#define gzopen(x,y) fopen(x,y)
#define gzgets(x,y,z) fgets(y,z,x)
#define gzclose(x) fclose(x)
#endif

#include <vector>
#include <algorithm>
using namespace std;

#include "chain.h"

#define NLINE    256
//0         1         2         3         4         5         6         7
//01234567890123456789012345678901234567890123456789012345678901234567890
//ATOM      1  N   LEU A   1      -3.862   4.293  48.792  1.00 24.16      1TCC 235
#define PDB_CHAIN 21
#define PDB_RESIDUE 17
#define PDB_RESIDUE_SIZE 3
#define PDB_RESIDUE_NUMBER 22
#define PDB_RESIDUE_NUMBER_SIZE 4
#define PDB_ICODE 26
#define PDB_ATOM 13
#define PDB_TFACTOR 60
#define PDB_TFACTOR_SIZE 6
#define PDB_PTS 30
#define PDB_PTS_SIZE 8
#define PDB_RESOLUTION 23

Aminoacid Chain::aa[21];		//amino acid array

Chain::Chain( char symbol, char *file ) {
  n_aa = 21;	for(int ii = 0; ii < n_aa; ii++) aa[ii].Set(ii);

  n_residues = 0;
  short_residues = 0;
  n_atoms = 0;
  n_models = 0;
  pdb_file = file;
  pdb_code[0] = 0;
  char *pos;
  for( pos = pdb_file + strlen( pdb_file ); pos > pdb_file + 3; pos-- ) 
    if( *pos == '.' ) {
      strncpy( pdb_code, pos - 4, 4 ); pdb_code[4] = 0;
      break;
    }
  if( symbol == '_' ) symbol = ' ';
  chain_symbol = symbol;

  resolution = 99.99;
}

Chain::Chain( Chain *c ) {
  n_aa = 21;	for(int ii = 0; ii < n_aa; ii++) aa[ii].Set(ii);

  n_models = c->n_models;
  n_residues = c->n_residues;
  short_residues = c->short_residues;
  n_atoms = c->n_atoms;
  chain_symbol = c->chain_symbol;
  resolution = c->resolution;
  pdb_file = c->pdb_file;
  for( int i = 0; i < 5; i++ ) pdb_code[i] = c->pdb_code[i];

  sequence = c->sequence;
  residue_chain = c->residue_chain;
  residue_model = c->residue_model;
  residue_number = c->residue_number;
  residue_icode = c->residue_icode;
  ica = c->ica;
  iresidues = c->iresidues;
  for( int i = 0; i < iresidues.size(); i++ ) {
    int i_residue = iresidues[i];
    residue_name.push_back( aa[i_residue].AA3 );
  }
  for( int i = 0; i < 3; i++ ) pts[i] = c->pts[i];
  tfactor = c->tfactor;
  faked_atom = c->faked_atom;
}

Chain::~Chain(void) {
}


int Chain::ReadPDB(void) {
  int print = 0;
  print = 0;
  int i, i_residue, i_atom, i_find, k, res_number, i_sum = 0;
  char res_num[PDB_RESIDUE_NUMBER_SIZE+1];
  int atom_size, res_atom_size, current_atom = -1, current_residue = -1, loaded_atoms, current_model;
  char current_chain;
  for(i = 0; i < NATM_MAX; i++) find_atm[i] = false;
#ifdef __unix__
  gzFile inp;
#else
  FILE *inp;
#endif
  char line[NLINE], line1[NLINE], tmp[NLINE], pa[5], icode;

  line[0] = 0;
  inp = gzopen(pdb_file, "r");
  if(!inp) {return false;}

  n_residues = 0; n_atoms = 0;
  i_residue = -1; i_atom = 0;
  loaded_atoms = 0; current_model = 0;

  pa[0] = pa[2] = pa[3] = 0;
  while(gzgets(inp, line, NLINE )){                                    // loop over file lines
    if(!strncmp(line, "REMARK   2 RESOLUTION", 21))	{
      if( strncmp(line+PDB_RESOLUTION, "NOT", 3 )	)
	resolution = atof(line+PDB_RESOLUTION);
      else
	resolution = -1.;
      continue;
    }
    if(!strncmp(line, "ENDMDL", 6)) {
      n_models++;
      continue;
    }
    if( ( strlen( line ) == 4 && !strncmp(line, "END", 3)) || !strncmp(line, "END   ", 6))
      break;

    if(strncmp(line, "ATOM", 4) && !( !strncmp(line, "HETATM", 6) && !strncmp(line+PDB_RESIDUE, "MSE", PDB_RESIDUE_SIZE) ) ) continue;
    if(chain_symbol != '*') {
      if(line[PDB_CHAIN] != chain_symbol) continue;
    }
    if(line[PDB_ATOM-1] == 'H') continue;
    if(line[PDB_ATOM] == 'H') continue;

    strncpy(tmp, &line[PDB_RESIDUE_NUMBER], PDB_RESIDUE_NUMBER_SIZE);
    tmp[PDB_RESIDUE_NUMBER_SIZE] = 0;
    res_number = atoi(tmp);

//incomplete residue
    if(i_residue >= 0 && ( res_number != current_residue || icode != line[PDB_ICODE] || current_chain != line[PDB_CHAIN] )) {
      FinalizeResidue( i_atom, i_residue, loaded_atoms );
      strncpy(res_num, &line1[PDB_RESIDUE_NUMBER], PDB_RESIDUE_NUMBER_SIZE+1);
      short_list.push_back( current_residue );
      short_residues++;
      if(print > 0) 
	fprintf( stderr, "%s:\tshort:%s", pdb_code, line1 );
    } 

    current_residue = res_number;
    icode = line[PDB_ICODE];
    current_chain = line[PDB_CHAIN];

    if(n_residues > 0) {
      strncpy(tmp, &line[PDB_RESIDUE_NUMBER], PDB_RESIDUE_NUMBER_SIZE+1);
      tmp[PDB_RESIDUE_NUMBER_SIZE+1] = 0;
      if( !strncmp(res_num, tmp, PDB_RESIDUE_NUMBER_SIZE+1) && icode == line[PDB_ICODE])
	{if(print > 1) fprintf( stderr, "%s:\ttail: %s", pdb_code, line ); continue;} 
    }
//residue type
    if(i_residue < 0) {
      for(i = 0; i < n_aa; i++) {
	if(!strncmp(&line[PDB_RESIDUE], aa[i].AA3, PDB_RESIDUE_SIZE)) {i_residue = i; break;}
      }
    }
    if(i_residue < 0)
      {if(print > 0) fprintf( stderr, "%s:\tresid:%s", pdb_code, line ); continue;} //Bad pdb string

    i_find = -1;
    for(i = 0; i < aa[i_residue].size; i++) {
      atom_size = 3;
      if(line[PDB_ATOM + 2] == ' ') atom_size = 2;
      if(line[PDB_ATOM + 1] == ' ') atom_size = 1;
      res_atom_size = strlen(aa[i_residue].atoms[i]);
      if(atom_size != res_atom_size) continue;

      for(k = 0; k < atom_size; k++) {
	if(line[PDB_ATOM + k] != aa[i_residue].atoms[i][k]) break;
	if(k == atom_size-1) {i_find = i; a_idx.push_back(i); break;} //find atom
      }
      if(i_find >= 0) break;
    }
    if(i_find == -1)
      {if(print > 0) fprintf( stderr, "%s:\tatom: %s", pdb_code, line ); continue;} //Bad string
    if( i_find == 6 && i_residue == 20 ) i_residue = 12;  // trick with 20 -> 12 is the workaround for MSE

    current_atom = loaded_atoms + i_find;
    if(find_atm[i_find] == false) {
      for(i = 0; i < 3; i++) {
	strncpy(tmp, &line[PDB_PTS + i*PDB_PTS_SIZE], PDB_PTS_SIZE);
	tmp[PDB_PTS_SIZE] = 0;
	tmp_pts[i].push_back( float(atof(tmp)) );
      }
      if(strlen(line) > (PDB_TFACTOR + PDB_TFACTOR_SIZE)) {
	strncpy(tmp, &line[PDB_TFACTOR], PDB_TFACTOR_SIZE);
	tmp[PDB_TFACTOR_SIZE] = 0;
	tmp_tfactor.push_back( float(atof(tmp)) );
      }
      else
        tmp_tfactor.push_back( 0. );
 
      if(i_find == 0) {
	ica.push_back( current_atom ); //index of N atom (first atom of residue)
	current_model = n_models;
	strncpy(tmp, &line[PDB_RESIDUE_NUMBER], PDB_RESIDUE_NUMBER_SIZE);
	tmp[PDB_RESIDUE_NUMBER_SIZE] = 0;
	residue_number.push_back( atoi(tmp) );
	residue_icode.push_back( line[PDB_ICODE] );
	sequence.push_back( aa[i_residue].AA1[0] );
	residue_chain.push_back( line[PDB_CHAIN] );
	residue_model.push_back( current_model );
	residue_name.push_back( aa[i_residue].AA3 );
	if( i_residue == 20 )  
	  iresidues.push_back( 12 );
	else
	  iresidues.push_back( i_residue );
      }
    }
    else
      a_idx.pop_back();
    for(k = 0; k < 4; k++) pa[k] = line[PDB_ATOM + k];
    if(find_atm[i_find] == false) {find_atm[i_find] = true; n_atoms++; i_atom++;}
    if(i_atom >= aa[i_residue].size) {
// residue is read; let's move data from temporary to permanent storage
      FinalizeResidue( i_atom, i_residue, loaded_atoms );
      strncpy(res_num, &line[PDB_RESIDUE_NUMBER], PDB_RESIDUE_NUMBER_SIZE+1);
    }
    strcpy( line1, line );
  }
  gzclose(inp);

  i = 1;
  if( a_idx.size() > 0 ) { // Wery last incomplete residue maybe still not in the permanent storage
    if( residue_number.size() > 0 )
      short_list.push_back( residue_number.back() );
    else
      short_list.push_back( -1 );
    short_residues++;
    if(print > 0) 
      fprintf( stderr, "%s:\tshort:%s", pdb_code, line1 );
    i = FinalizeResidue( i_atom, i_residue, loaded_atoms );
  }

  if( n_models == 0 ) n_models = 1;
  n_atoms = pts[0].size();
  return true;
}

int Chain::FinalizeResidue( int &i_atom, int &i_residue, int &loaded_atoms ) {
  int i, k;

  if( !find_atm[0] || !find_atm[1] ) { // N atom is used for "ica". We cannot continue unless it has been defined. We also require CA atom to exist
    if( find_atm[0] ) { // we cannot proceed without CA, but if N atom was found, we need to pop residue-related data
      ica.pop_back();
      residue_number.pop_back();
      residue_icode.pop_back();
      sequence.pop_back();
      residue_chain.pop_back();
      residue_model.pop_back();
      residue_name.pop_back();
      iresidues.pop_back();
    }
    n_atoms -= a_idx.size();
    for(i = 0; i < NATM_MAX; i++) find_atm[i] = false;
    tmp_pts[0].clear();  tmp_pts[1].clear();  tmp_pts[2].clear();
    tmp_tfactor.clear();
    a_idx.clear();
    i_residue = -1;
    i_atom = 0;
  }
  else {
// working around MODELLER atom order
    for( i = 0; i < aa[i_residue].size; i++ ) {
      tfactor.push_back(999.);
      faked_atom.push_back(true);
      for( k = 0; k < 3; k++ )
	pts[k].push_back(0.);
    }
    for( i = 0; i < a_idx.size(); i++ ) {
      if( find_atm[a_idx[i]] ) {
	tfactor[loaded_atoms+a_idx[i]] = tmp_tfactor[i];
	faked_atom[loaded_atoms+a_idx[i]] = not find_atm[a_idx[i]];
	for( k = 0; k < 3; k++ )
	  pts[k][loaded_atoms+a_idx[i]] = tmp_pts[k][i];
      }
    }
    tmp_tfactor.clear();
    for( k = 0; k < 3; k++ )
      tmp_pts[k].clear();
    a_idx.clear();

    loaded_atoms += aa[i_residue].size;
    n_residues ++; i_residue = -1; i_atom = 0;
    for(i = 0; i < NATM_MAX; i++) find_atm[i] = false;
  }
  return true;
}


int Chain::WritePDB(char *pdb_file) {
  int i, j, ia, na, i_seq = 0, i_res; char tmp[255], mdl, chain_smb;
  FILE *out;
  out = fopen(pdb_file, "wt");
  if(!out) {return false;}
  mdl = -1;

  for(i = 0; i < n_residues; i++) {
    if( n_models > 1 )
      if( mdl != residue_model[i] ) {
	mdl = residue_model[i];
        fprintf( out, "MODEL      %3d\n", mdl+1 );
      }
    ia = ica[i]; i_res = iresidues[i]; na = aa[i_res].size;
    chain_smb = residue_chain[i];
    for(j = 0; j < na; j++) {
      i_seq ++;
      if( pts[0][ia+j] == 0. && pts[1][ia+j] == 0. && pts[2][ia+j] == 0. && tfactor[ia+j] == 999. )
	continue;
      sprintf(tmp, "ATOM %6d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f  1.00%6.2f",
	      i_seq, aa[i_res].atoms[j], aa[i_res].AA3, chain_smb, residue_number[i],
	      pts[0][ia+j], pts[1][ia+j], pts[2][ia+j], tfactor[ia+j]);
      fprintf(out, "%s\n", tmp);
    }
    if( n_models > 1 )
      if( i+1 < n_residues )
        if( mdl != residue_model[i+1] )
	  fprintf( out, "ENDMDL\n" );
  }
//0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 
//ATOM    157  CD2 LEU A  19       1.640  27.834  22.650  1.00 13.64           C  
  if( n_models > 1 )
    fprintf( out, "ENDMDL\n" );
  fclose(out);
  return true;
}

Residue** Chain::MkRes(void) {
  int i, i1, j, j1;
  Residue **r = new Residue*[n_residues];

  for( i = 0; i < n_residues; i++ )
    r[i] = new Residue( this, i );
  for( i = 1; i < n_residues; i++ )
    if( r[i-1]->chain_name == r[i]->chain_name && r[i-1]->model == r[i]->model ) {
      r[i]->prev = r[i-1];
      r[i-1]->next = r[i];
    } else
      r[i]->prev = r[i-1]->next = NULL;
  r[0]->prev = r[n_residues-1]->next = NULL;

  return r;
}
