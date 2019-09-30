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
#include "../forcefield.h"

int main(int argc, char **argv)
{
  if(argc >= 2){

    ForceField ff;
    LoadFFParams(argv[argc-1], &ff);

    MOLECULE molecule;
    for(size_t i = 1; i < argc-1; i++){
      NewMOL2Molecule(&molecule, argv[i]);
      AtomAnalyzer(&molecule, 5);
      printf("Processing: %s\n", molecule.molname);
      //PrintMoleculeAtomFingerprint(molecule);
      FFParamsAssigner(molecule, &ff);
      DelMolecule(&molecule);
    }

    SaveFFParams(ff, argv[argc-1]);
    DeleteForceField(&ff);
  }
  else{
    printf("\nUsage %s file.mol2 file2.mol2 ... fileN.mol2 ffield.txt\n\n", argv[0]);
    return -1;
  }

  return 0;
}
