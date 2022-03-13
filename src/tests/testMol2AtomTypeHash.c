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
  printf("Hybridization types: \n");
  printf("%d %d %d %d %d %d %d %d %d\n", S, SP, SP2, SP3, SP3D, SP3D2, AR, OTHER, UNSPECIFIED);
  for(int i = 0; i < molecule.n_atoms; i++){
    printf("id: %d ", i+1);
    int t = molecule.atoms[i].ainfo.hybrid;
    switch(t){
      case 1:
        printf("S ");
        break;
      case 2:
        printf("SP ");
        break;
      case 3:
        printf("SP2 ");
        break;
      case 4:
        printf("SP3 ");
        break;
      case 5:
        printf("SP3D ");
        break;
      case 6:
        printf("SP3D2 ");
        break;
      case 7:
        printf("AR ");
        break;
      case 8:
        printf("OTHER ");
        break;
      default:
        printf("UNSPECIFIED ");
    }

    printf("fp: %s\n",molecule.atoms[i].ainfo.atype_hash);
  }
  printf("END\n");
  DelMolecule(&molecule);
  return 0;
}
