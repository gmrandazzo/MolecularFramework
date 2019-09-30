/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../molecule.h"
#include "../mol3Daligner.h"

/* Test to check if hash are correctly defined! */

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2 cx cy cz filetransl.mol2\n\n", argv[0]);
    return -1;
  }
  printf("%s\n", argv[1]);
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  TranslateConformation(molecule, atof(argv[2]), atof(argv[3]), atof(argv[4]));
  SaveMol2Molecule(molecule,  argv[5]);
  PrintMolecule(molecule);
  DelMolecule(&molecule);
  return 0;
}
