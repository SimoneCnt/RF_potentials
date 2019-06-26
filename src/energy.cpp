#include <stdio.h>
#include <string.h>

#include "residue.h"
#include "energy.h"

#ifndef __unix__
#include <limits>
#define INFINITY numeric_limits<double>::infinity()
#endif

#include "point_3.h"
//#include "version.cpp"


Energy::Energy( char *PotenFile ) {
// bool
  Directed = Processed = Smooth = false;
// double
  aa_ttl = pair_total = BinSize = MinDist = AverDist = MaxDist = AppCutoff = E = 0.;
// int
  Nprot = AverType = n_atoms = Nres = NatmType = natm = ChainSepMin= ChainSepMax = Nbins =  ChainDependent =  DepthDep = 0;

  pdb = NULL;
  r = NULL;

  poten = NULL;
  energy = NULL;

  Descr[0] = 0;



  na[0] = na[4] = 20; na[1] =167; na[2] = na[3] = 0;

  bool state;
  for( int read_count = 0; read_count < 3; read_count++ ) {
    state = ReadPotential( PotenFile );
    if( state )
      break;
    sleep( 5 );
  }
  if( !state ) {
    printf( "can not read \"%s\"!\n", PotenFile );
    exit( 1 );
  }

  BinSize = (MaxDist - MinDist) / Nbins;
}

double Energy::Execute( Residue **R ) {
  int i, j, res00, res0, res11, res1;

  r = R;

  pdb = r[0]->chain;
  n_atoms = pdb->pts[0].size();
  Nres = pdb->n_residues;

  bool GlyCbNeeded = false;
  if( NatmType == 4 ) GlyCbNeeded = true;
  if(  GlyCbNeeded )
    for( i = 0; i < pdb->n_residues; i++ )
      r[i]->MakeVirtualCb();


  energy = new double*[Nres];
  for( i = 0; i < Nres; i++ ) {
    energy[i] = new double[Nbins];
    for( j = 0; j < Nbins; j++ )
      energy[i][j] = 0.;
  }
//                  A  R   N   D   C   Q   E   G   H   I   L   K   M   F    P    S    T    W    Y    V
  int Rshft[20] = { 0, 5, 16, 24, 32, 38, 47, 56, 60, 70, 78, 86, 95, 103, 114, 121, 127, 134, 148, 160 };
  int natm_ref[MAXATMTYPES][NATM];
  for( i = 0; i < MAXATMTYPES; i++ )
    for( j = 0; j < NATM; j++ )
      natm_ref[i][j] = -1;


  for( j = 0; j < 20; j++ ) {
    natm_ref[0][Rshft[j]+1] = j;	// CA
    if( j == 7 )
      natm_ref[4][Rshft[j]+1] = j;	// CA for Gly
    else
      natm_ref[4][Rshft[j]+4] = j;	// CB
  }
  for( j = 0; j < NATM; j++ )
    natm_ref[1][j] = j;			// all Heavy Atoms

  int CSepMax = ( ChainSepMax > Nres ) ? Nres : ChainSepMax;
  Atom *a0, *a1;

  if( ChainSepMax == -1 )
    CSepMax = Nres;

  for( i = 0; i < Nres; i++ ) {
    Point_3 a;
    if( Directed ) {    
      a.x = (r[i]->atoms[4]->x[0] - r[i]->atoms[1]->x[0]);
      a.y = (r[i]->atoms[4]->x[1] - r[i]->atoms[1]->x[1]);
      a.z = (r[i]->atoms[4]->x[2] - r[i]->atoms[1]->x[2]);
      a/= ~a;
    }
    for( int i1 = 0; i1 < r[i]->n_atoms; i1++ ) {
      res00 = natm_ref[NatmType][Rshft[r[i]->res_idx]+i1];
      if( res00 == -1 ) continue;
      if( NatmType == 4 && r[i]->res_idx == 7 ) // for CB potentials we use artificial CB now
	a0 = r[i]->atoms[4];
      else
	a0 = r[i]->atoms[i1];
      if( a0 == NULL ) continue;
      if( NatmType != 4 || r[i]->res_idx != 7 ) // for CB potentials we use artificial CB now
	if( a0->faked ) continue;

      int j0 = (ChainDependent%2) ? i+ChainSepMin : i-CSepMax;	//  odd means assymmetric
      for( j = j0; j <= i + CSepMax; j++ ) {
	if( j < 0 || j >= Nres ) continue;
	int k = i-j;
	if( k < ChainSepMin && k > -ChainSepMin )
	  continue;

	for( int j1 = 0; j1 < r[j]->n_atoms; j1++ ) {
	  res11 = natm_ref[NatmType][Rshft[r[j]->res_idx]+j1];
	  if( res11 == -1 ) continue;
	  if( NatmType == 4 && r[j]->res_idx == 7 ) // for CB potentials we use artificial CB now
	    a1 = r[j]->atoms[4];
	  else
	    a1 = r[j]->atoms[j1];
	  if( a1 == NULL )
	    continue;
	  if( NatmType != 4 || r[j]->res_idx != 7 ) // for CB potentials we use artificial CB now
	    if( a1->faked ) continue;

	  if( r[i]->chain_name != r[j]->chain_name && (int)(ChainDependent/2) )    // ChainDependent > 1 means "chain only"
	    continue;

	  if( Directed ) {
	    Point_3 b( (r[j]->atoms[4]->x[0] - r[j]->atoms[1]->x[0]), (r[j]->atoms[4]->x[1] - r[j]->atoms[1]->x[1]), (r[j]->atoms[4]->x[2] - r[j]->atoms[1]->x[2]) ); b /= ~b;
	    Point_3 c( (r[j]->atoms[1]->x[0] - r[i]->atoms[1]->x[0]), (r[j]->atoms[1]->x[1] - r[i]->atoms[1]->x[1]), (r[j]->atoms[1]->x[2] - r[i]->atoms[1]->x[2]) ); c /= ~c;
	    if( a*b < 0 )
	      if( a*c > 0 ) {
		res0 = res00 + na[NatmType];
		res1 = res11 + na[NatmType];
	      }
	      else {
		res0 = res00 + 2*na[NatmType];
		res1 = res11 + 2*na[NatmType];
	      }
	    else {
	      res0 = res00;
	      res1 = res11;
	    }
	  }
	  else {
	    res0 = res00;
	    res1 = res11;
	  }


	  double dd = a1->Distance( a0 );
	  double ddd = dd - MinDist;
	  int d = (int)( ddd / BinSize );

	  if( d < 0 ) d = 0;          // this is bad: if i==j we score pseudo-contact to itself!!!
	  if( d >= Nbins ) continue;

	  if( poten[res0][res1][d] != INFINITY )
	    energy[i][d] += poten[res0][res1][d];
	  else energy[i][d] += MY_INF;
	}
      }
    }
  }
  E = SumUp( 0., AppCutoff );
  return E;
}


#include <sys/stat.h>
#include <time.h>
bool Energy::ReadPotential( const char *FileName ) {  // read binary potential
  int i, j, k, ver;
  char buff[1024];
  FILE *f;

  int VERSION = 17;

  Directed = false;

  struct stat filestat;
  time_t t20070101 = 1167627600; // Mon Jan  1 00:00:00 2007
  stat( FileName, &filestat );
  bool OldFile = ( filestat.st_mtime < t20070101 ) ? true : false;

  int EVP_SAVED_MD_SIZE = 36; // EVP_MAX_MD_SIZE used before Jan 2007
  if( EVP_SAVED_MD_SIZE > EVP_MAX_MD_SIZE ) {
    printf( "OpenSSL conflict: EVP_SAVED_MD_SIZE = %d > EVP_MAX_MD_SIZE = %d !\n", EVP_SAVED_MD_SIZE, EVP_MAX_MD_SIZE );
    exit( 1 );
  }
  f = fopen( FileName, "rb" );
  if( !f ) {
    printf( "PairPotential \"%s\" not found!\n", FileName );
    exit( 1 );
  }
  fread_unlocked( Descr, sizeof(char), 128, f );
  fread_unlocked( &ver, sizeof(int), 1, f );
  if( ver != VERSION ) {
    printf( "\007\n  Error! The potential file \"%s\" has obsolete format version:\n    \"%d\" found, whereas current one is \"%d\"\n  Please update the \"%s\" file and run the program again.\n", FileName, ver, VERSION, FileName );
    exit( 1 );
  }
  fread_unlocked( &Nprot, sizeof(int), 1, f );
  fread_unlocked( &NatmType, sizeof(int), 1, f );
  if( NatmType < 0 ) {
    Directed = true;
    NatmType = -NatmType;
  }
  if( NatmType >= MAXATMTYPES ) {
    printf( "AtomType = %d, > than MaxAtmtypes = %d\n", NatmType, MAXATMTYPES );
    exit( 1 );
  }
  natm = Directed ? 3*na[NatmType] : na[NatmType];
  fread_unlocked( &ChainSepMin, sizeof(int), 1, f );
  fread_unlocked( &ChainSepMax, sizeof(int), 1, f );
  fread_unlocked( &ChainDependent, sizeof(int), 1, f );
  fread_unlocked( &MinDist, sizeof(double), 1, f );
  fread_unlocked( &AverDist, sizeof(double), 1, f );
  fread_unlocked( &MaxDist, sizeof(double), 1, f );
  fread_unlocked( &Nbins, sizeof(int), 1, f );
  fread_unlocked( &AverType, sizeof(int), 1, f );
  fread_unlocked( &AppCutoff, sizeof(double), 1, f );
  fread_unlocked( &DepthDep, sizeof(int), 1, f );
  fread_unlocked( &Smooth, sizeof(bool), 1, f );
  fread_unlocked( &Processed, sizeof(bool), 1, f );
  if( OldFile )
    fread_unlocked( md_value, sizeof(unsigned char), EVP_SAVED_MD_SIZE, f );

  poten = new double**[natm];
  for( i = 0; i < natm; i++ ) {
    poten[i] = new double*[natm];
    for( j = 0; j < natm; j++ ) {
      poten[i][j] = new double[Nbins];
      fread_unlocked( poten[i][j], sizeof(double), Nbins, f );
    }
  }
  fread_unlocked( &aa_ttl, sizeof(double), 1, f );
  fread_unlocked( &pair_total, sizeof(double), 1, f );
  if( !OldFile ) {
    fread_unlocked( &EVP_SAVED_MD_SIZE, sizeof(int), 1, f );
    if( EVP_SAVED_MD_SIZE > EVP_MAX_MD_SIZE ) {
      printf( "OpenSSL conflict: EVP_SAVED_MD_SIZE = %d > EVP_MAX_MD_SIZE = %d !\n", EVP_SAVED_MD_SIZE, EVP_MAX_MD_SIZE );
      exit( 1 );
    }
    fread_unlocked( md_value, sizeof(unsigned char), EVP_SAVED_MD_SIZE, f );
  }

  if( ferror( f ) ) {
    fclose( f );
    printf( "! ! ! ! ! ! !       \"ferror\" error detected while reading \"%s\"!\n",  FileName );
    return( false );
  } else
    fclose( f );

  return( checksum_verify() );
}


bool Energy::checksum_verify(void) {
  bool is_OK = true;
  unsigned char my_md_value[EVP_MAX_MD_SIZE];
  unsigned int md_len = checksum(my_md_value);

  for( unsigned int i = 0; i < md_len; i++ ) 
    if( md_value[i] != my_md_value[i] ) {
      is_OK = false;
      unsigned int ii;
      printf( "# PairPotential reading error!   SHA1 signatures for the \"poten\" are different!\n#      Stored one is " );
      PrintSignature();
      printf( "\n#      Calcuated is  " );
      PrintSignature( stdout, my_md_value );
      printf( "\n#\n" );
      break;
    }
  return( is_OK );
}


unsigned int Energy::checksum(unsigned char *md_val) {
  EVP_MD_CTX *mdctx = EVP_MD_CTX_create();
  const EVP_MD *md;
  unsigned int md_len;

  OpenSSL_add_all_digests();
//  md = EVP_get_digestbyname("sha1");
  md = EVP_sha1();

//  EVP_MD_CTX_init(&mdctx);
//  EVP_DigestInit_ex(&mdctx, md, NULL);
  EVP_DigestInit(mdctx, md);
  for( int i = 0; i < natm; i++ )
    for( int j = 0; j < natm; j++ )
      EVP_DigestUpdate(mdctx, poten[i][j], Nbins*sizeof(double));
  EVP_DigestFinal( mdctx, md_val, &md_len );
//  EVP_DigestFinal_ex( &mdctx, md_val, &md_len );
  EVP_MD_CTX_cleanup(mdctx);
EVP_MD_CTX_destroy(mdctx);
  return( md_len );
}

void Energy::PrintSignature( FILE *out ) { PrintSignature( out, md_value ); }
void Energy::PrintSignature( FILE *out, unsigned char *md_val ) {
  for( unsigned int i = 0; i < 20; i++ ) fprintf( out,"%02x", md_val[i]);
}

void Energy::PrintPotential( FILE *out, int mode ) {
  int i, j, k;
  int natm1 = Directed ? natm / 3 : natm;


    const char *Aa1[MAXATMTYPES][NATM] = { 
      { "ALA.CA", "ARG.CA", "ASN.CA", "ASP.CA", "CYS.CA", "GLN.CA", "GLU.CA", "GLY.CA", "HIS.CA", "ILE.CA", "LEU.CA", "LYS.CA", "MET.CA", "PHE.CA", "PRO.CA", "SER.CA", "THR.CA", "TRP.CA", "TYR.CA", "VAL.CA" },
      { "ALA.N", "ALA.CA", "ALA.C", "ALA.O", "ALA.CB", "ARG.N", "ARG.CA", "ARG.C", "ARG.O", "ARG.CB", "ARG.CG", "ARG.CD", "ARG.NE", "ARG.CZ", "ARG.NH1", "ARG.NH2", "ASN.N", "ASN.CA", "ASN.C", "ASN.O", "ASN.CB", "ASN.CG", "ASN.OD1", "ASN.ND2", "ASP.N", "ASP.CA", "ASP.C", "ASP.O", "ASP.CB", "ASP.CG", "ASP.OD1", "ASP.OD2", "CYS.N", "CYS.CA", "CYS.C", "CYS.O", "CYS.CB", "CYS.SG", "GLN.N", "GLN.CA", "GLN.C", "GLN.O", "GLN.CB", "GLN.CG", "GLN.CD", "GLN.OE1", "GLN.NE2", "GLU.N", "GLU.CA", "GLU.C", "GLU.O", "GLU.CB", "GLU.CG", "GLU.CD", "GLU.OE1", "GLU.OE2", "GLY.N", "GLY.CA", "GLY.C", "GLY.O", "HIS.N", "HIS.CA", "HIS.C", "HIS.O", "HIS.CB", "HIS.CG", "HIS.ND1", "HIS.CD2", "HIS.CE1", "HIS.NE2", "ILE.N", "ILE.CA", "ILE.C", "ILE.O", "ILE.CB", "ILE.CG1", "ILE.CG2", "ILE.CD1", "LEU.N", "LEU.CA", "LEU.C", "LEU.O", "LEU.CB", "LEU.CG", "LEU.CD1", "LEU.CD2", "LYS.N", "LYS.CA", "LYS.C", "LYS.O", "LYS.CB", "LYS.CG", "LYS.CD", "LYS.CE", "LYS.NZ", "MET.N", "MET.CA", "MET.C", "MET.O", "MET.CB", "MET.CG", "MET.SD", "MET.CE", "PHE.N", "PHE.CA", "PHE.C", "PHE.O", "PHE.CB", "PHE.CG", "PHE.CD1", "PHE.CD2", "PHE.CE1", "PHE.CE2", "PHE.CZ", "PRO.N", "PRO.CA", "PRO.C", "PRO.O", "PRO.CB", "PRO.CG", "PRO.CD", "SER.N", "SER.CA", "SER.C", "SER.O", "SER.CB", "SER.OG", "THR.N", "THR.CA", "THR.C", "THR.O", "THR.CB", "THR.OG1", "THR.CG2", "TRP.N", "TRP.CA", "TRP.C", "TRP.O", "TRP.CB", "TRP.CG", "TRP.CD1", "TRP.CD2", "TRP.NE1", "TRP.CE2", "TRP.CE3", "TRP.CZ2", "TRP.CZ3", "TRP.CH2", "TYR.N", "TYR.CA", "TYR.C", "TYR.O", "TYR.CB", "TYR.CG", "TYR.CD1", "TYR.CD2", "TYR.CE1", "TYR.CE2", "TYR.CZ", "TYR.OH", "VAL.N", "VAL.CA", "VAL.C", "VAL.O", "VAL.CB", "VAL.CG1", "VAL.CG2" },
      { "unused" },
      { "unused" },
      { "ALA.CB", "ARG.CB", "ASN.CB", "ASP.CB", "CYS.CB", "GLN.CB", "GLU.CB", "GLY.CA", "HIS.CB", "ILE.CB", "LEU.CB", "LYS.CB", "MET.CB", "PHE.CB", "PRO.CB", "SER.CB", "THR.CB", "TRP.CB", "TYR.CB", "VAL.CB" }
    };

  const char *AverTypeDescript[2] = { "AVER_PSEUDOCOUNT", "AVER_SIPPLSCHEME" };  // list of possible averaging types
  const char *AtomTypeDescript[2*MAXATMTYPES] = { "CA", "AllAtoms", "unused", "unused", "CB",
				       "unused", "unused", "unused", "unused", "CB_Directed" };

  fprintf( out, "#\n# PairPotential \"%s\" of \"%s\" atom type:\n", Descr, AtomTypeDescript[(int)Directed*MAXATMTYPES + NatmType] ); fflush( stdout );
  fprintf( out, "# Chain Separation Min = %d\n", ChainSepMin ); fflush( stdout );
  fprintf( out, "# Chain Separation Max = %d\n#\n", ChainSepMax ); fflush( stdout );
  fprintf( out, "# Minimal Distance = %4.1f\n", MinDist ); fflush( stdout );
  fprintf( out, "# Maximal Distance = %4.1f\n", MaxDist ); fflush( stdout );
  fprintf( out, "# Number of bins = %d\n#\n", Nbins ); fflush( stdout );
  fprintf( out, "# Averaging type =\t %s\n", AverTypeDescript[AverType] );
  fprintf( out, "# Maximal distance of the potential applicability = %4.1f\n", AppCutoff );
  fprintf( out, "# ChainDependent = %d\n", ChainDependent );
  
  fprintf( out, "# Built using %d proteins\n", Nprot ); fflush( stdout );
  fprintf( out, "# Total number of the atoms %0.0f\n", aa_ttl );
  fprintf( out, "# Total number of the all atom pairs  %0.0f\n#\n", pair_total );
  fprintf( out, "# SHA1 signature for the \"poten\" is " );
  PrintSignature( out );
  fprintf( out, "\n#\n" );
/* */
  if( mode > 0 ) {
    for( k = 0; k < Nbins; k++ ) {
      fprintf( out, "\nBin %d (%4.1f A < x < %4.1f A)\n", k, ( (float)k * (MaxDist - MinDist) / Nbins + MinDist ),
	       ( (float)(k+1) * (MaxDist - MinDist) / Nbins + MinDist ) );
      if( Directed )
	fprintf( out, " j parallel i\n" );
      fflush( out );
      for( j = 0; j < natm1; j++ ) {
        fprintf( out, "\t%s", Aa1[NatmType][j] ); fflush( out );
      }
      fprintf( out, "\n" ); fflush( out );
      for( i = 0; i < natm1; i++ ) {
        fprintf( out, "%s\t", Aa1[NatmType][i] ); fflush( out );
        for( j = 0; j < natm1; j++ ) {
          fprintf( out, "%8.3f", poten[i][j][k] ); fflush( out );
	}
        fprintf( out, "\n" ); fflush( out );
      }


      if( Directed ) {
	fprintf( out, " j anti-parallel i, facing each other\n" );
	for( j = 0; j < natm1; j++ ) {
	  fprintf( out, "\t%s", Aa1[NatmType][j] ); fflush( out );
	}
	fprintf( out, "\n" ); fflush( out );
	for( i = natm1; i < 2*natm1; i++ ) {
	  fprintf( out, "%s\t", Aa1[NatmType][i-natm1] ); fflush( out );
	  for( j = natm1; j < 2*natm1; j++ ) {
	    fprintf( out, "%8.3f", poten[i][j][k] ); fflush( out );
	  }
	  fprintf( out, "\n" ); fflush( out );
	}

	fprintf( out, " j anti-parallel i, looking opposite directions\n" );
	for( j = 0; j < natm1; j++ ) {
	  fprintf( out, "\t%s", Aa1[NatmType][j] ); fflush( out );
	}
	fprintf( out, "\n" ); fflush( out );
	for( i = 2*natm1; i < natm; i++ ) {
	  fprintf( out, "%s\t", Aa1[NatmType][i-2*natm1] ); fflush( out );
	  for( j = 2*natm1; j < natm; j++ ) {
	    fprintf( out, "%8.3f", poten[i][j][k] ); fflush( out );
	  }
	  fprintf( out, "\n" ); fflush( out );
	}
      }
    }
  }
}



double Energy::SumUp( double DistMin, double DistMax ) {
  int binmin, binmax;
  binmin = BinMin( DistMin );
  binmax = BinMax( DistMax );
  if( binmin < 0 ) binmin = 0;
  if( binmax > Nbins ) binmax = Nbins;

  double e = 0.;

  for( int d = binmin; d < binmax; d++ ) 
    for( int i = 0; i < Nres; i++ )
      if( energy[i][d] != INFINITY )
	e += energy[i][d];
      else e += MY_INF;
  return e;
}

Energy::Energy( const Energy &e ) {
  int i, j;

  pdb = e.pdb;
  n_atoms = e.n_atoms;
  Nbins = e.Nbins;
  r = e.r;
  MinDist = e.MinDist;
  MaxDist = e.MaxDist;
  BinSize = e.BinSize;
  Nprot = e.Nprot;
  Processed = e.Processed;
  for( i = 0; i < EVP_MAX_MD_SIZE; i++ )
    md_value[i] = e.md_value[i];
  aa_ttl = e.aa_ttl;
  pair_total = e.pair_total;
  AverType = e.AverType;
  Nres = e.Nres;
  NatmType = e.NatmType;
  Directed = e.Directed;
  natm = e.natm;
  ChainSepMin = e.ChainSepMin;
  ChainSepMax = e.ChainSepMax;
  AverDist = e.AverDist;
  AppCutoff = e.AppCutoff;
  ChainDependent = e.ChainDependent;
  DepthDep = e.DepthDep;
  Smooth = e.Smooth;
  E = e.E;

  energy = new double*[Nres];
  for( i = 0; i < Nres; i++ ) {
    energy[i] = new double[Nbins];
    for( j = 0; j < Nbins; j++ )
      energy[i][j] = e.energy[i][j];
  }

  poten = new double**[natm];
  for( i = 0; i < natm; i++ ) {
    poten[i] = new double*[natm];
    for( j = 0; j < natm; j++ ) {
      poten[i][j] = new double[Nbins];
      for( int k = 0; k < Nbins; k++ )
	poten[i][j][k] = e.poten[i][j][k];
    }
  }
}

Energy::~Energy(void) {
  if( energy != NULL ) {
    for( int i = 0; i < Nres; i++ ) {
      delete [] energy[i];
      energy[i] = NULL;
    }
    delete [] energy;
    energy = NULL;
  }

  if( poten != NULL ) {
    for( int i = 0; i < natm; i++ ) {
      for( int j = 0; j < natm; j++ ) {
	if( poten[i][j] != NULL ) {
	  delete [] poten[i][j];
	  poten[i][j] = NULL;
	}
      }
      delete [] poten[i];
      poten[i] = NULL;
    }
    delete [] poten;
    poten = NULL;
  }
}
