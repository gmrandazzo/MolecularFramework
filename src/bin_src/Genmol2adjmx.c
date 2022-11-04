/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Dfields.h"



int main(int argc, char **argv)
{
  if(argc != 3){
    printf("\nUsage %s file.mol2 [0: bond length adj, 1: bond color adj]\n\n", argv[0]);
    return -1;
  }

  size_t i, j;
  MOLECULE molecule;
  // printf("Load molecule %s ...\n", molecule.molname);
  NewMOL2Molecule(&molecule, argv[1]);
  matrix *adjmx;
  initMatrix(&adjmx);
  if(atoi(argv[2]) == 0){
    BondLenghtAdjMatGen(&molecule, adjmx);
  }
  else{
    BondColorAdjMatGen(&molecule, adjmx);
  }

  printf("atom types,");
  for(i = 0; i < molecule.n_atoms-1; i++){
    printf("%s,", molecule.atoms[i].type);
  }
  printf("%s\n", molecule.atoms[molecule.n_atoms-1].type);

  for(i = 0; i < molecule.n_atoms; i++){
    printf("%s,", molecule.atoms[i].type);
    for(j = 0; j < molecule.n_atoms-1; j++){
      printf("%.4f,", adjmx->data[i][j]);
    }
    printf("%.4f\n", adjmx->data[i][molecule.n_atoms-1]);
  }

  DelMatrix(&adjmx);
  DelMolecule(&molecule);
  return 0;
}
