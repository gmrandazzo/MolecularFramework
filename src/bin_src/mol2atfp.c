/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include "../molecule.h"
#include "../atomanalysis.h"


int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2 bond_depth\n\n", argv[0]);
    return -1;
  }
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, atoi(argv[2]));
  /* Header
   * Molname_atomname_hybrid charge C N O S P H F Cl Br I S SP SP2 SP3 SP3D SP3D2 AR OTHER >>
   * >> N * (c_sp c_sp2 c_sp3 c_ar n_sp n_sp2 n_ar o_sp2 o_sp3 o_ar s_sp2 s_sp3 s_sp3d s_sp3d2 s_ar p f cl br i other)
   */

  for(int i = 0; i < molecule.n_atoms; i++){
    printf("%s_%d_%s_", molecule.molname, i+1, molecule.atoms[i].asymbl);
    //printf("%s ", molecule.atoms[i].type);

    int t = molecule.atoms[i].ainfo.hybrid;
    switch(t){
      case 1:
        printf("S,");
        break;
      case 2:
        printf("SP,");
        break;
      case 3:
        printf("SP2,");
        break;
      case 4:
        printf("SP3,");
        break;
      case 5:
        printf("SP3D,");
        break;
      case 6:
        printf("SP3D2,");
        break;
      case 7:
        printf("AR,");
        break;
      case 8:
        printf("OTHER,");
        break;
      default:
        printf("UNSPECIFIED,");
    }

    printf("%.4f,",molecule.atoms[i].charge);

    //"C N O S P H F Cl Br I"
    if(strcmp(molecule.atoms[i].asymbl, "C") == 0){
      printf("1,0,0,0,0,0,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "N") == 0){
      printf("0,1,0,0,0,0,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "O") == 0){
      printf("0,0,1,0,0,0,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "S") == 0){
      printf("0,0,0,1,0,0,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "P") == 0){
      printf("0,0,0,0,1,0,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "H") == 0){
      printf("0,0,0,0,0,1,0,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "F") == 0){
      printf("0,0,0,0,0,0,1,0,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "Cl") == 0){
      printf("0,0,0,0,0,0,0,1,0,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "Br") == 0){
      printf("0,0,0,0,0,0,0,0,1,0,");
    }
    else if(strcmp(molecule.atoms[i].asymbl, "I") == 0){
      printf("0,0,0,0,0,0,0,0,0,1,");
    }
    else{
      printf("99,99,99,99,99,99,99,99,99,99 ");
    }

    /* S, SP, SP2, SP3, SP3D, SP3D2, OTHER, UNSPECIFIED   */
    if(molecule.atoms[i].ainfo.hybrid == 1){
      printf("1,0,0,0,0,0,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 2){
      printf("0,1,0,0,0,0,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 3){
      printf("0,0,1,0,0,0,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 4){
      printf("0,0,0,1,0,0,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 5){
      printf("0,0,0,0,1,0,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 6){
      printf("0,0,0,0,0,1,0,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 7){
      printf("0,0,0,0,0,0,1,0,0,");
    }
    else if(molecule.atoms[i].ainfo.hybrid == 8){
      printf("0,0,0,0,0,0,0,1,0,");
    }
    else{
      printf("0,0,0,0,0,0,0,0,1");
    }


    /*int c = 0;
    while(molecule.atoms[i].ainfo.atype_hash[c] != '\0'){
      //printf("%c ", molecule.atoms[i].ainfo.atype_hash[c]);
      c++;
    }
    printf("%d \n", c);*/
    printf("%s\n", molecule.atoms[i].ainfo.atype_hash);

  }
  DelMolecule(&molecule);
  return 0;
}
