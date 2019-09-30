/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../cppwrap.h"
#include "../miscdesc.h"

/* Test to check if hash are correctly defined! */

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2\n\n", argv[0]);
    return -1;
  }
  printf("%s\n", argv[1]);
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  //PrintMolecule(molecule);
  PrintMoleculeRings(molecule);
  printf("Formal Charge %d\n", GetFormalCharge(molecule));
  //SaveMol2Molecule(molecule, "suca.mol2");
  DelMolecule(&molecule);
  return 0;
}
