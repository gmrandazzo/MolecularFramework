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
#include "../miscdesc.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.sdf\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewSDFMolecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  PrintMolecule(molecule);
  DelMolecule(&molecule);
  return 0;
}
