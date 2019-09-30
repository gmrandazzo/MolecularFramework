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
#include "../atomdesc.h"
#include "../periodic_table.h"

/* Test to check if hash are correctly defined! */

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2 file.pqr\n\n", argv[0]);
    return -1;
  }
  //printf("%s\n", argv[1]);
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);

  /*
  ForceField ff;
  LoadFFParams(argv[2], &ff);
  AssignParams(&molecule, ff, vanderwaals);
  */
  for(int i = 0; i < molecule.n_atoms; i++){
    molecule.atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule.atoms[i].asymbl);
  }

  printf("%s Sum of charges: %.8f\n", molecule.molname, GetTotalCharge(molecule));

  SavePQRFile(molecule, argv[2]);
  //DeleteForceField(&ff);
  DelMolecule(&molecule);
  return 0;
}
