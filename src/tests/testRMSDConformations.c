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
#include "../mol3Daligner.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf(" RMSD between two 3D conformations A to a molecule B\n");
    printf("\nUsage %s confA.mol2 confB.mol2\n\n", argv[0]);
    return -1;
  }

  MOLECULE A, B;
  NewMOL2Molecule(&A, argv[1]);
  NewMOL2Molecule(&B, argv[2]);

  double rmsd = ComputeRMSD(A, B);
  printf("RMSD: %.3f\n", rmsd);

  DelMolecule(&A);
  DelMolecule(&B);
  return 0;
}
