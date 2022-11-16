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
    printf(" Align a 3D molecule A to a molecule B\n");
    printf("\nUsage %s align_A.mol2 in_B.mol2\n\n", argv[0]);
    return -1;
  }

  MOLECULE A, B;
  NewMOL2Molecule(&A, argv[1]);
  AtomAnalyzer(&A, 1);
  NewMOL2Molecule(&B, argv[2]);
  AtomAnalyzer(&B, 1);
  double rmsd = Align3DConformations(A, B);
  //double rmsd = SimplexAlign3DConformations(A, B);
  printf("RMSD: %.3f\n", rmsd);

  SaveMol2Molecule(A,  argv[1]);
  DelMolecule(&A);
  DelMolecule(&B);
  return 0;
}
