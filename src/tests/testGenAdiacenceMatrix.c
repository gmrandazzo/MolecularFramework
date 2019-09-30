/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <scientific.h>
#include "../molecule.h"
#include "../atomanalysis.h"


int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2\n\n", argv[0]);
    return -1;
  }
  printf("%s\n", argv[1]);
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  matrix *adjmx;
  initMatrix(&adjmx);
  AdjMatGen(&molecule, &adjmx);
  PrintMatrix(adjmx);
  DelMatrix(&adjmx);
  DelMolecule(&molecule);
  return 0;
}
