/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../molecule.h"
#include "../cppwrap.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2 ID \n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  int *atom_depth;
  getAtomDepth(&molecule, atoi(argv[2]), 3, &atom_depth);
  for(int j = 0; j < molecule.n_atoms; j++){
    printf("%d %d\n", j+1, atom_depth[j]);
  }
  DelMolecule(&molecule);
  free(atom_depth);
  return 0;
}
