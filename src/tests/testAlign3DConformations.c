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
  MOLECULE A, B;
  NewMOL2Molecule(&A, argv[1]);
  AtomAnalyzer(&A, 1);
  NewMOL2Molecule(&B, argv[2]);
  AtomAnalyzer(&B, 1);
  double rmsd = Align3DConformations(A, B);

  printf("TEST ALIGN A 3D CONFORMER A TO B: ");
  if(FLOAT_EQ(rmsd, 0.f, 1E-2)){
    printf("PASS\n");
  }
  else{
    printf("ERROR!\n");
    abort();
  }

  DelMolecule(&A);
  DelMolecule(&B);
  return 0;
}
