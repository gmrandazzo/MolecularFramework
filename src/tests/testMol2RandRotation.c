/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../molecule.h"
#include "../mol3Daligner.h"

int main(int argc, char **argv)
{
  if(argc == 3){
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    RandomConformationRotation(&molecule);
    SaveMol2Molecule(molecule, argv[2]);
    DelMolecule(&molecule);
    return 0;
  }
  else{
    printf("\nUsage %s infile.mol2 outfile.mol2 \n\n", argv[0]);
    return -1;
  }
}
