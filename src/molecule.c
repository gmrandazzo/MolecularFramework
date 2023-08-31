/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "version.h"
#include "molecule.h"
#include "misc.h"
#include "atomanalysis.h"


struct enumtypes
{
   char *str;
};

struct enumtypes mol2field[] = {
  {"@<TRIPOS>ALT_TYPE"},
  {"@<TRIPOS>ANCHOR_ATOM"},
  {"@<TRIPOS>ASSOCIATED_ANNOTATION"},
  {"@<TRIPOS>CENTER_OF_MASS"},
  {"@<TRIPOS>CENTROID"},
  {"@<TRIPOS>COMMENT"},
  {"@<TRIPOS>CRYSIN"},
  {"@<TRIPOS>DICT"},
  {"@<TRIPOS>DATA_FILE"},
  {"@<TRIPOS>EXTENSION_POINT"},
  {"@<TRIPOS>FF_PBC"},
  {"@<TRIPOS>FFCON_ANGLE"},
  {"@<TRIPOS>FFCON_DIST"},
  {"@<TRIPOS>FFCON_MULTI"},
  {"@<TRIPOS>FFCON_RANGE"},
  {"@<TRIPOS>FFCON_TORSION"},
  {"@<TRIPOS>LINE"},
  {"@<TRIPOS>LSPLANE"},
  {"@<TRIPOS>NORMAL"},
  {"@<TRIPOS>QSAR_ALIGN_RULE"},
  {"@<TRIPOS>RENDERING_ATTRS"},
  {"@<TRIPOS>RING_CLOSURE"},
  {"@<TRIPOS>ROTATABLE_BOND"},
  {"@<TRIPOS>SEARCH_DIST"},
  {"@<TRIPOS>SEARCH_OPTS"},
  {"@<TRIPOS>SET"},
  {"@<TRIPOS>SUBSTRUCTURE"},
  {"@<TRIPOS>U_FEAT"},
  {"@<TRIPOS>UNITY_ATOM_ATTR"},
  {"@<TRIPOS>UNITY_BOND_ATTR"}
};

static int skipMol2Field(char *str){
  int i;
  const int sz = sizeof(mol2field) / sizeof(mol2field[0]);

  for(i = 0; i < sz; i++){
    if(strstr(str, mol2field[i].str) != NULL){
      return 0;
    }
  }
  return -1;
}

static inline void Mol2NameFill(int i, int j, char *tl, MOLECULE *molecule)
{
  if(i == 0 && j == 0){
    strcpy(molecule->molname, tl);
    // sscanf(tl, "%s", molecule->molname);
    //      printf("Molecular Name: %s\n",tl);
  }
  else if(i == 1 && j == 0){
  //     printf("Num of atoms: %s\n",tl);
    molecule->n_atoms = atoi(tl);
  }
  else if( i == 1 && j == 1){
//     printf("Num of bondss: %s\n",tl);
    molecule->n_bonds = atoi(tl);
  }
}

void getAtomSymbol(char *atomtype, char **atomsymbl)
{
  (*atomsymbl) = strdup(atomtype);
  (*atomsymbl) = strtok((*atomsymbl), ".");
}

static inline void Mol2AtomFill( int i, int j, char *tl, MOLECULE *molecule)
{
  if ( j == 1 ){
//     printf("Atom Name: %s\n",tl);
    strcpy(molecule->atoms[i].name, tl);
    /*sscanf(tl, "%s", molecule->atoms[i].name);*/
  }
  else if( j == 2 ){
//     printf("x : %s\n",tl);
    molecule->atoms[i].coord.x = atof(tl);
  }
  else if( j == 3 ){
//     printf("y : %s\n",tl);
    molecule->atoms[i].coord.y = atof(tl);
  }
  else if( j == 4 ){
//     printf("y : %s\n",tl);
    molecule->atoms[i].coord.z = atof(tl);
  }
  else if( j == 5 ){
//     printf("atom type/symbol : %s\n",tl);
    strcpy(molecule->atoms[i].type, tl);
    char *atomsymbl = NULL;
    getAtomSymbol(tl, &atomsymbl);
    strcpy(molecule->atoms[i].asymbl, atomsymbl);
    free(atomsymbl);
    /*sscanf(tl, "%s", molecule->atoms[i].type);*/
  }
  else if( j == 6 ){
    molecule->atoms[i].subst_id = atoi(tl);
  }
  else if( j == 7 ){
    strcpy(molecule->atoms[i].subst_name, tl);
  }
  else if( j == 8 ){
    molecule->atoms[i].charge = atof(tl);
   }
   else{
     return;
   }
}

static inline void Mol2BondFill( int i, int j, char *tl, MOLECULE *molecule )
{
  if ( j == 1 ){
//     printf("bond origin id: %s\n",tl);
    molecule->bonds[i].origin_atom_id = atoi(tl)-1;
  }
  else if( j == 2){
//     printf("bond target id: %s\n",tl);
    molecule->bonds[i].target_atom_id = atoi(tl)-1;
  }
  else if( j == 3 ){
//     printf("bond type: %s\n",tl);
    /*sscanf(tl, "%s", molecule->bonds[i].type);*/
    if(strcmp(tl, "am") == 0){
      strcpy(molecule->bonds[i].type, "1");
    }
    else{
      strcpy(molecule->bonds[i].type, tl);
    }

  }
  else{
    return;
  }
}

static void ReadMOL2(const char *filename, MOLECULE *molecule)
{
  int i, j;
  char l[LINE_MAX], *lPtr, *tl;
  FILE *fmol;

  int /*molname = 1,*/ atom = 1, bond = 1;

 if((fmol = fopen (filename, "r")) == NULL){
    printf ("File %s not found!!!", filename);
    return;
  }

  i = 0;
  while(fgets(l, LINE_MAX, fmol) != NULL){
    if(strstr(l, "#") != NULL){
      continue;
    }
    /*else if(strstr(l, "@<TRIPOS>MOLECULE") != NULL){
      i = 0;
      molname = 0;
      atom = 1;
      bond = 1;
      continue;
    }*/
    else if(strstr(l, "@<TRIPOS>ATOM") != NULL){
      i = 0;
      //molname = 1;
      atom = 0;
      bond = 1;
      continue;
    }
    else if(strstr(l, "@<TRIPOS>BOND") != NULL){
      i = 0;
      //molname = 1;
      atom = 1;
      bond = 0;
      continue;
    }
    else if(skipMol2Field(l) == 0){
      //molname = 1;
      atom = 1;
      bond = 1;
      continue;
    }
    else{
      /*Tokenize
      for(tl = strtok(Trim(l), " \t"), j = 0; tl; tl = strtok (NULL, " \t"), j++){
        if(molname == 0){
          //Mol2NameFill( i, j, tl, molecule);
          continue;
        }
        else if(atom == 0){
          Mol2AtomFill(i, j, tl, molecule);
        }
        else if(bond == 0){
          Mol2BondFill( i, j, tl, molecule);
        }
        else
          continue;
      }*/
      //char s[2] = " ";
      lPtr = l;
      Trim(lPtr);
      if(atom == 0 && bond == 1){
        j = 0;
        while((tl = strsep(&lPtr, " ")) != NULL ){
          if(strcmp(tl, "") == 0){
            continue;
          }
          else{
            //printf("atom %d <%s>\n", j, tl);
            Mol2AtomFill(i, j, tl, molecule);
            j++;
          }
        }
      }
      else if(bond == 0 && atom == 1){
        j = 0;
        while((tl = strsep(&lPtr, " \t")) != NULL ){
          if(strcmp(tl, "") == 0){
            continue;
          }
          else{
            //printf("atom %d <%s>\n", j, tl);
            Mol2BondFill(i, j, tl, molecule);
            j++;
          }
        }
      }
      i++;
    }
  }
  fclose(fmol);
}

void InitializeMolecule(MOLECULE *molecule)
{
  for(size_t i = 0; i < molecule->n_atoms; i++){
    strcpy(molecule->atoms[i].name, "None");
    strcpy(molecule->atoms[i].type, "None");
    strcpy(molecule->atoms[i].asymbl, "None");
    strcpy(molecule->atoms[i].subst_name, "None");
    molecule->atoms[i].subst_id = -9999;
    molecule->atoms[i].coord.x = -9999; molecule->atoms[i].coord.y = -9999; molecule->atoms[i].coord.z = -9999;
    molecule->atoms[i].charge = -9999;
    molecule->atoms[i].radius = -9999;
    molecule->atoms[i].ainfo.atype_hash = NULL;
    molecule->atoms[i].ainfo.valence_electrons = -9999;
    molecule->atoms[i].ainfo.lonepairs = -9999;
    molecule->atoms[i].ainfo.n_Others = -9999;
    molecule->atoms[i].ainfo.n_Cl = -9999;
    molecule->atoms[i].ainfo.n_F = -9999;
    molecule->atoms[i].ainfo.n_Br = -9999;
    molecule->atoms[i].ainfo.n_I = -9999;
    molecule->atoms[i].ainfo.n_P = -9999;
    molecule->atoms[i].ainfo.n_S = -9999;
    molecule->atoms[i].ainfo.n_N = -9999;
    molecule->atoms[i].ainfo.n_O = -9999;
    molecule->atoms[i].ainfo.n_C = -9999;
    molecule->atoms[i].ainfo.n_H = -9999;
    molecule->atoms[i].ainfo.aromatic = -9999;
    molecule->atoms[i].ainfo.hybrid = UNSPECIFIED;
    molecule->atoms[i].ainfo.nsp = -9999;
    molecule->atoms[i].ainfo.nsp2 = -9999;
    molecule->atoms[i].ainfo.nsp3 = -9999;
    molecule->atoms[i].ainfo.nar = -9999;
    molecule->atoms[i].ainfo.stereo = -9999;
    molecule->atoms[i].ainfo.cycle = -9999;
    molecule->atoms[i].ainfo.connectivity = -9999;
  }

  for(size_t i = 0; i < molecule->n_bonds; i++){
    strcpy(molecule->bonds[i].type, "None");
    molecule->bonds[i].origin_atom_id = -9999;
    molecule->bonds[i].target_atom_id = -9999;
  }

  molecule->rings = NULL;
}

void NewEmptyMolecule(MOLECULE *molecule , int n_atoms, int n_bonds)
{
  molecule->n_atoms = n_atoms;
  molecule->n_bonds = n_bonds;
  if(n_atoms > 0)
    molecule->atoms = malloc(sizeof(ATOM)*molecule->n_atoms);
  if(n_bonds > 0)
    molecule->bonds = malloc(sizeof(BOND)*molecule->n_bonds);
  InitializeMolecule(molecule);
}


void NewMOL2Molecule(MOLECULE *molecule , const char *filename)
{
  int i, j;
  char l[LINE_MAX], *lPtr, *tl;
  FILE *fmol;
  int molname = 1;

  if((fmol = fopen (filename, "r")) == NULL){
    printf ("File %s not found!!!", filename);
    return;
  }

  i = 0;
  while(fgets(l, LINE_MAX, fmol) != NULL){
    if (strstr(l, "#") != NULL) { // skip line
      continue;
    }
    else if(strstr(l, "@<TRIPOS>MOLECULE") != NULL){ // skip this line and set to fill the next line value after her
      i = 0;
      molname = 0;
      continue;
    }
    else if(strstr(l, "@<TRIPOS>ATOM") != NULL){ // skip this line and set to fill the next line value after her
      molname = 1;
      continue;
    }
    else if(strstr(l, "@<TRIPOS>BOND") != NULL){ // skip this line and set to fill the next line value after her
      i = 0;
      molname = 1;
      continue;
    }
    else if(skipMol2Field(l) == 0){
      molname = 1;
      continue;
    }
    else{
      /* Tokenize
      for (tl = strtok(Trim(l), " \t"), j = 0; tl; tl = strtok (NULL, " \t"), j++){
        if(molname == 0){
          Mol2NameFill(i, j, tl, molecule);
        }
        else
          continue;
      }*/
      if(molname == 0){
        lPtr = l;
        Trim(lPtr);
        if(i == 0){
          //strcpy(molecule->molname, lPtr);
          //sscanf(lPtr, "%s", molecule->molname);
          Mol2NameFill(0, 0, lPtr, molecule);
        }
        else{
          j = 0;
          while((tl = strsep(&lPtr, " ")) != NULL ){
            if(strcmp(tl, "") == 0){
              continue;
            }
            else{
              //printf("atom %d <%s>\n", j, tl);
              Mol2NameFill(i, j, tl, molecule);
              j++;
            }
          }
        }
      }
      i++;
    }
  }

  molecule->atoms = malloc(sizeof(ATOM)*molecule->n_atoms);
  molecule->bonds = malloc(sizeof(BOND)*molecule->n_bonds);
  /*Initialize values*/
  InitializeMolecule(molecule);

  fclose(fmol);
  ReadMOL2(filename, molecule);

  /*Atom analysis to assign atom type
    AtomAnalyzer(molecule, 1);
  */
}

void SaveMol2Molecule(MOLECULE molecule, const char *filename)
{
  int i;
  FILE *fmol;

  if((fmol = fopen (filename, "w")) == NULL){
    printf ("Unable to save file %s\n", filename);
    return;
  }
  fprintf(fmol, "@<TRIPOS>MOLECULE\n");
  fprintf(fmol, "%s\n", molecule.molname);
  //" 13 12 0 0 0"
  //fprintf(fmol, "%3d %3d 0 0 0\n", molecule.n_atoms, molecule.n_bonds);
  fprintf(fmol, " %d %d 0 0 0\n", molecule.n_atoms, molecule.n_bonds);
  fprintf(fmol, "SMALL\n");
  fprintf(fmol, "GASTEIGER\n\n@<TRIPOS>ATOM\n");
/*
# Example of format ATOM
#      1 CB         -0.6670   -1.1470   -0.3838 C.3     1  ALA1       -0.0395
#     11 HXT         1.9717    0.1248   -1.1260 H       1  ALA1        0.2951
*/
  for(i = 0; i < molecule.n_atoms; i++){
    fprintf(fmol, "%7d %-4s       % 3.4f   % 3.4f   % 3.4f %-4s %4d  %4s       % 2.4f\n",
            i+1,
            molecule.atoms[i].name,
            molecule.atoms[i].coord.x,
            molecule.atoms[i].coord.y,
            molecule.atoms[i].coord.z,
            molecule.atoms[i].type,
            molecule.atoms[i].subst_id,
            molecule.atoms[i].subst_name,
            molecule.atoms[i].charge);
  }
  fprintf(fmol, "@<TRIPOS>BOND\n");

/* Example of format BOND
#     1     2     1    1
#    41    21    29   ar
*/
  for(i = 0; i < molecule.n_bonds; i++){
    fprintf(fmol, "%6d %5d %5d %4s\n",
            i+1,
            molecule.bonds[i].origin_atom_id+1,
            molecule.bonds[i].target_atom_id+1,
            molecule.bonds[i].type);
  }
  fclose(fmol);
}

void NewSDFMolecule(MOLECULE *molecule, const char *filesdf)
{
  size_t j = 0;
  char l[LINE_MAX], *lPtr, *tl;
  FILE *fmol;

  if((fmol = fopen (filesdf, "r")) == NULL){
    printf ("File %s not found!!!", filesdf);
    return;
  }
  int cc, ia, ib; // global counter and internal file counter for get atom field and bond field
  cc = ia = ib = 0;

  fseek (fmol, 0, SEEK_SET);
  while(fgets(l, LINE_MAX, fmol) != NULL){
    lPtr = l;
    Trim(lPtr);
    if(cc == 0){ // first line is the name
      if(strlen(lPtr) == 0){ // skip line
        char name[] = "No_named_molecule";
        copystr(name, molecule->molname);
      }
      else{
        copystr(lPtr, molecule->molname);
      }
//       printf("%s\n", molecule->molname);
    }
    else if(cc == 3){ // third line are number of atoms and number of bonds
      // get  natoms and n_bonds
      j = 0;
      for(tl = strtok(lPtr, " "); tl; tl = strtok (NULL, " ")){
        if(j == 0){
           molecule->n_atoms = atoi(tl);
        }
        else if(j == 1){
          molecule->n_bonds = atoi(tl);
        }
        else{
          break;
        }
        j++;
      }
//       printf("natoms %d nbonds %d\n", molecule->n_atoms, molecule->n_bonds);
      molecule->atoms = malloc(sizeof(ATOM)*molecule->n_atoms);
      molecule->bonds = malloc(sizeof(BOND)*molecule->n_bonds);
    }
    else if(cc>3 && ia < molecule->n_atoms){
      j = 0;
      for(tl = strtok(lPtr, " "); tl; tl = strtok (NULL, " ")){
        if(j == 0){
           molecule->atoms[ia].coord.x = atof(tl);
        }
        else if(j == 1){
          molecule->atoms[ia].coord.y = atof(tl);
        }
        else if(j == 2){
          molecule->atoms[ia].coord.z = atof(tl);
        }
        else if(j == 3){
          Trim(tl);
          copystr(tl, molecule->atoms[ia].type);
          copystr(tl, molecule->atoms[ia].asymbl);
          copystr(tl, molecule->atoms[ia].name);
        }
        else if(j == 6){
          int stereo = atoi(tl);
          if(stereo == 0)
            molecule->atoms[ia].ainfo.stereo = 0;
          else if(stereo == 1) // odd
            molecule->atoms[ia].ainfo.stereo = 5;
          else if(stereo == 2) // even
            molecule->atoms[ia].ainfo.stereo = 6;
          else
            molecule->atoms[ia].ainfo.stereo = 7; // undefined

          molecule->atoms[ia].charge = 0.f;
        }
        else{
          /*notthing*/
        }
        j++;
      }

//       printf("%s %f %f %f %d\n", molecule->atoms[ia].type, molecule->atoms[ia].coord.x, molecule->atoms[ia].coord.y, molecule->atoms[ia].coord.z, molecule->atoms[ia].stereo);
      ia++;
    }
    else if(cc > 3 && ia >=molecule->n_atoms && ib<molecule->n_bonds){
      j = 0;
      for(tl = strtok(lPtr, " "); tl; tl = strtok (NULL, " ")){
        if(j == 0){
           molecule->bonds[ib].origin_atom_id = atoi(tl)-1;
        }
        else if(j == 1){
          molecule->bonds[ib].target_atom_id = atoi(tl)-1;
        }
        else if(j == 2){
          Trim(tl);
          copystr(tl, molecule->bonds[ib].type);
        }
        else{
          break;
        }
        j++;
      }

//       printf("%d %d %s\n", molecule->bonds[ib].origin_atom_id, molecule->bonds[ib].target_atom_id, molecule->bonds[ib].type);
      ib++;
    }
    cc++;
  }
  fclose(fmol);

  /*Atom analysis to assign atom type
  AtomAnalyzer(molecule, 1);
  */
}

void SaveSDFMolecule(MOLECULE *molecule, const char *filesdf){ printf("Function not yet implemented\n"); }

void DelMolecule(MOLECULE *molecule)
{
  for(size_t i = 0; i < molecule->n_atoms; i++){
    free(molecule->atoms[i].ainfo.atype_hash);
  }
  if(molecule->n_atoms > 0)
    free(molecule->atoms);
  if(molecule->n_bonds > 0)
    free(molecule->bonds);

  if(molecule->rings != NULL){
    for(int i = 0; i < molecule->n_rings; i++){
      free(molecule->rings[i].atoms);
    }
    free(molecule->rings);
  }
}


void SavePQRFile(MOLECULE molecule, const char *filename)
{
  int i;
  FILE *fmol;

  if((fmol = fopen (filename, "w")) == NULL){
    printf ("Unable to save file %s\n", filename);
    return;
  }
  int major, minor, patch;
  get_version(&major, &minor, &patch);

  fprintf(fmol, "REMARK   1 PQR file generated by MolDesc version %d %d %d\n", major, minor, patch);
  fprintf(fmol, "REMARK   1\n");
  fprintf(fmol, "REMARK   1 MOLECULE: %s\n", molecule.molname);
  fprintf(fmol, "REMARK   1\n");
/*
# Example of format ATOM
#ATOM      3  H   UNL A   1      -5.193   1.648  -0.714  0.2249 1.2000
#...
*/
  for(i = 0; i < molecule.n_atoms; i++){
    fprintf(fmol, "ATOM %6d %2s   UNL A   1      % 5.3f  % 5.3f  % 5.3f % 5.4f %5.4f\n",
            i+1,
            molecule.atoms[i].name,
            molecule.atoms[i].coord.x,
            molecule.atoms[i].coord.y,
            molecule.atoms[i].coord.z,
            molecule.atoms[i].charge,
            molecule.atoms[i].radius);
  }
  fprintf(fmol, "TER\nEND\n");
  fclose(fmol);
}

void PrintMolecule(MOLECULE molecule)
{
  int i;
  printf("molname: %s\n", molecule.molname);
  printf("N°. Atoms: %d\n", molecule.n_atoms);
  printf("N°. Bonds: %d\n", molecule.n_bonds);

  puts("Atoms");
  for(i = 0 ; i < molecule.n_atoms; i++){
    char hybrd[4];
    if(molecule.atoms[i].ainfo.hybrid == SP){
      strcpy(hybrd, "SP");
    }
    else if(molecule.atoms[i].ainfo.hybrid == SP2){
      strcpy(hybrd, "SP2");
    }
    else if(molecule.atoms[i].ainfo.hybrid == SP3){
      strcpy(hybrd, "SP3");
    }
    //else if(molecule.atoms[i].ainfo.hybrid == AR){
    else{
      strcpy(hybrd, "Ar");
    }

    printf("%d\n", i+1);
    printf("atom: %s\n", molecule.atoms[i].name);
    printf("file_atomtype: %s\n", molecule.atoms[i].type);
    printf("atom_symbol: %s\n", molecule.atoms[i].asymbl);
    printf("x: %f y: %f z: %f\n", molecule.atoms[i].coord.x, molecule.atoms[i].coord.y, molecule.atoms[i].coord.z);
    printf("partial_charge: %f\n", molecule.atoms[i].charge);
    printf("lonepairs: %d\n", molecule.atoms[i].ainfo.lonepairs);
    printf("hybridisation: %s\n", hybrd);
    printf("stereochemistry: %d\n", molecule.atoms[i].ainfo.stereo);
    printf("iscycle: %d\n", molecule.atoms[i].ainfo.cycle);
    printf("isaromatic: %d\n", molecule.atoms[i].ainfo.aromatic);
    printf("n. SP: %d\n", molecule.atoms[i].ainfo.nsp);
    printf("n. SP2: %d\n", molecule.atoms[i].ainfo.nsp2);
    printf("n. SP3: %d\n", molecule.atoms[i].ainfo.nsp3);
    printf("n. AR: %d\n", molecule.atoms[i].ainfo.nar);
    printf("n. H: %d\n", molecule.atoms[i].ainfo.n_H);
    printf("n. C: %d\n", molecule.atoms[i].ainfo.n_C);
    printf("n. O: %d\n", molecule.atoms[i].ainfo.n_O);
    printf("n. N: %d\n", molecule.atoms[i].ainfo.n_N);
    printf("n. S: %d\n", molecule.atoms[i].ainfo.n_S);
    printf("n. P: %d\n", molecule.atoms[i].ainfo.n_P);
    printf("n. F: %d\n", molecule.atoms[i].ainfo.n_F);
    printf("n. Cl: %d\n", molecule.atoms[i].ainfo.n_Cl);
    printf("n. Br: %d\n", molecule.atoms[i].ainfo.n_Br);
    printf("n. I: %d\n", molecule.atoms[i].ainfo.n_I);
    printf("n. Others: %d\n", molecule.atoms[i].ainfo.n_Others);
    printf("-------------------------------------------------------\n");
  }

  puts("Bonds");
  for(i = 0; i < molecule.n_bonds; i++){
    printf("%d %d %d %s\n", i+1, molecule.bonds[i].origin_atom_id, molecule.bonds[i].target_atom_id, molecule.bonds[i].type);
  }
}

void PrintMoleculeRings(MOLECULE molecule)
{
  size_t i, j;
  printf("\nNumber of rings found: %d\n", molecule.n_rings);
  for(i = 0 ; i < molecule.n_rings; i++){
    for(j = 0; j < molecule.rings[i].size-1; j++){
      printf("%d ", molecule.rings[i].atoms[j]+1);
    }
    printf("%d\n", molecule.rings[i].atoms[molecule.rings[i].size-1]+1);
  }
}


void PrintMoleculeAtomInfo(MOLECULE molecule)
{
  int i;
  printf("molname: %s\n", molecule.molname);
  printf("N°. Atoms: %d\n", molecule.n_atoms);
  printf("N°. Bonds: %d\n", molecule.n_bonds);

  puts("Atoms");
  puts("ID   Name  X  Y  Z  Type  Charge  Stereo  NHydrogen  Hybrid  Cycle  nsp  nsp2  nsp3  nar  connectivity");
  for(i = 0 ; i < molecule.n_atoms; i++){
    printf("%3d %5s %5.3lf %5.3lf %5.3lf %5s %3.3lf %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n",
          i+1,
          molecule.atoms[i].name,
          molecule.atoms[i].coord.x,
          molecule.atoms[i].coord.y,
          molecule.atoms[i].coord.z,
          molecule.atoms[i].type,
          molecule.atoms[i].charge,
          molecule.atoms[i].ainfo.stereo,
          molecule.atoms[i].ainfo.n_H,
          molecule.atoms[i].ainfo.hybrid,
          molecule.atoms[i].ainfo.cycle,
          molecule.atoms[i].ainfo.lonepairs,
          molecule.atoms[i].ainfo.nsp,
          molecule.atoms[i].ainfo.nsp2,
          molecule.atoms[i].ainfo.nsp3,
          molecule.atoms[i].ainfo.nar,
          molecule.atoms[i].ainfo.connectivity);
  }

  puts("Bonds");
  for(i = 0; i < molecule.n_bonds; i++){
    printf("%d %d %d %s\n", i+1, molecule.bonds[i].origin_atom_id, molecule.bonds[i].target_atom_id, molecule.bonds[i].type);
  }
}

void PrintMoleculeAtomFingerprint(MOLECULE molecule)
{
  int i;
  for(i = 0 ; i < molecule.n_atoms; i++){
    printf("%s %s %.4f %.4f %.4f %.4f Atom %d in molecule %s\n", molecule.atoms[i].ainfo.atype_hash, molecule.atoms[i].type, 0.f, 0.f, 0.f, 0.f, i+1, molecule.molname);
  }
}
