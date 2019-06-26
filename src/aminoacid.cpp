/*
"G","GLY", 4,"Glycine"		N CA C O
"A","ALA", 5,"Alanine"		N CA C O CB
"C","CYS", 6,"Cysteine"		N CA C O CB SG
"S","SER", 6,"Serine"		N CA C O CB OG
"T","THR", 7,"Threonine"	N CA C O CB OG1 CG2
"P","PRO", 7,"Proline"		N CA C O CB CG  CD
"V","VAL", 7,"Valine"		N CA C O CB CG1 CG2
"L","LEU", 8,"Leucine"		N CA C O CB CG  CD1 CD2
"I","ILE", 8,"Isoleucine"	N CA C O CB CG1 CG2 CD1
"N","ASN", 8,"Asparagine"	N CA C O CB CG  OD1 ND2
"D","ASP", 8,"Aspartic acid"	N CA C O CB CG  OD1 OD2
"M","MET", 8,"Methionine",	N CA C O CB CG  SD  CE
"E","GLU", 9,"Glutamic acid"	N CA C O CB CG  CD  OE1 OE2
"Q","GLN", 9,"Glutamine"	N CA C O CB CG  CD  OE1 NE2
"K","LYS", 9,"Lysine"		N CA C O CB CG  CD  CE  NZ
"H","HIS",10,"Histidine"	N CA C O CB CG  ND1 CD2 CE1 NE2
"F","PHE",11,"Phenylalanine"	N CA C O CB CG  CD1 CD2 CE1 CE2 CZ
"R","ARG",11,"Arginine"		N CA C O CB CG  CD  NE  CZ  NH1 NH2
"Y","TYR",12,"Tyrosine"		N CA C O CB CG  CD1 CD2 CE1 CE2 CZ  OH
"W","TRP",14,"Tryptophan"	N CA C O CB CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2
"M","MSE", 8,"SelenoMethionine",N CA C O CB CG SE   CE
*/

#include <string.h>
#include "aminoacid.h"

char AA1List[NAA+1][2] = {  "A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",
			    "L",  "K",  "M",  "F",  "P",  "S",  "T",  "W",  "Y",  "V", "M"};
char AA3List[NAA+1][4] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
			  "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL", "MSE"};
int sizeList[NAA+1]	 = {   5,    11,   8,    8,    6,    9,    9,    4,    10,   8,
			       8,     9,   8,   11,    7,    6,    7,   14,    12,   7, 8};

char atomsList[NAA+1][NATM_MAX][4] = {
	"N", "CA", "C", "O", "CB",   "",     "",    "",    "",    "",    "",    "",    "",    "",//ALA
	"N", "CA", "C", "O", "CB", "CG",   "CD",  "NE",  "CZ", "NH1", "NH2",	"",    "",    "",//ARG
	"N", "CA", "C", "O", "CB", "CG",  "OD1", "ND2",	   "",    "",    "",    "",    "",    "",//ASN
	"N", "CA", "C", "O", "CB", "CG",  "OD1", "OD2",	   "",    "",    "",    "",    "",    "",//ASP
	"N", "CA", "C", "O", "CB", "SG",     "",    "",    "",    "",    "",    "",    "",    "",//CYS
	"N", "CA", "C", "O", "CB", "CG",   "CD", "OE1", "NE2",    "",    "",    "",    "",    "",//GLN
	"N", "CA", "C", "O", "CB", "CG",   "CD", "OE1", "OE2",	  "",    "",    "",    "",    "",//GLU
	"N", "CA", "C", "O",   "",   "",     "",    "",    "",    "",    "",    "",    "",    "",//GLY
	"N", "CA", "C", "O", "CB", "CG",  "ND1", "CD2", "CE1", "NE2",    "",    "",    "",    "",//HIS
	"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1",	   "",    "",    "",    "",    "",    "",//ILE
	"N", "CA", "C", "O", "CB", "CG",  "CD1", "CD2",	   "",    "",    "",    "",    "",    "",//LEU
	"N", "CA", "C", "O", "CB", "CG",   "CD",  "CE",  "NZ",    "",    "",    "",    "",    "",//LYS
	"N", "CA", "C", "O", "CB", "CG",   "SD",  "CE",    "",    "",    "",    "",    "",    "",//MET
	"N", "CA", "C", "O", "CB", "CG",  "CD1", "CD2", "CE1", "CE2",  "CZ",    "",    "",    "",//PHE
	"N", "CA", "C", "O", "CB", "CG",   "CD",    "",    "",    "",    "",    "",    "",    "",//PRO
	"N", "CA", "C", "O", "CB", "OG",     "",    "",    "",    "",    "",    "",    "",    "",//SER
	"N", "CA", "C", "O", "CB", "OG1", "CG2",    "",	   "",    "",    "",    "",    "",    "",//THR
	"N", "CA", "C", "O", "CB", "CG",  "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",//TRP
	"N", "CA", "C", "O", "CB", "CG",  "CD1", "CD2", "CE1", "CE2",  "CZ",  "OH",    "",    "",//TYR
	"N", "CA", "C", "O", "CB", "CG1", "CG2",    "",    "",    "",    "",    "",    "",    "",//VAL
	"N", "CA", "C", "O", "CB", "CG",    "E",  "CE",    "",    "",    "",    "",    "",    "" //MSE
};


void Aminoacid::Set(int i) {
  int j, k;
  strncpy(AA1, AA1List[i], 2);
  strncpy(AA3, AA3List[i], 4);
  size = sizeList[i];
  for(j = 0; j < NATM_MAX; j++) strncpy(atoms[j], atomsList[i][j], 4);
}
