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
#include "../massanalysis.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2\n\n", argv[0]);
    return -1;
  }

  char *bform;
  double  exbfmw;
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  GetMolecularFormula(molecule, &bform);
  GetExactMWfromMolecularFormula(bform, &exbfmw);
  printf("%s;%.10f;0.0;0.0;%s\n", molecule.molname, exbfmw, bform);
  DelMolecule(&molecule);
  free(bform);
  return 0;
}
