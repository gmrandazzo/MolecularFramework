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
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  printf("%s %f\n", molecule.molname, GetTotalPartialCharge(molecule));
  //printf("Molecule: %s Charge: %d\n", molecule.molname, GetFormalCharge(molecule));
  //printf("Molecule: %s Charge: %d\n", molecule.molname, GetNetCharge(molecule));
  DelMolecule(&molecule);
  return 0;
}
