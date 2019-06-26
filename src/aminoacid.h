#ifndef _AMINOACID_H_
#define _AMINOACID_H_

#define NAA 20	//number of aminoacids
#define NBB 12	//max number of backbones
#define NATM_MAX 14
//#define N_GLY 7  //No CB atom

//char AA3List[NAA][4] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
//					    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
enum RESIDUE_NAME {N_ALA, N_ARG, N_ASN, N_ASP, N_CYS, N_GLN, N_GLU, N_GLY, N_HIS, N_ILE,
                   N_LEU, N_LYS, N_MET, N_PHE, N_PRO, N_SER, N_THR, N_TRP, N_TYR, N_VAL};
enum ATOM_NAME {NA_N, NA_CA, NA_C, NA_O, NA_CB};

class Aminoacid {
 public:
  char AA1[2];
  char AA3[4];
  int size;

  char atoms[NATM_MAX][4];

  Aminoacid(void) {};
  ~Aminoacid(void) {};
  void Set(int i);
};
#endif
