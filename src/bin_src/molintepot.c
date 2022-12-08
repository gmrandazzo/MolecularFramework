/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Dfields.h"

int main(int argc, char **argv)
{
  if(argc != 6){
    printf("\nUsage %s file.mol2 [charge] [N voxel points (INT)] [Angstrom grid size (INT)] outfield.dx\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  printf("Load molecule %s ...\n", molecule.molname);

  int grid_resolution = atoi(argv[3]);
  int grid_size = atoi(argv[4]);

  printf("Calculate Electrostatic Potential ...\n");
  VOXEL *field;
  VoxelElectrostaticPotentialInteractionCalculator(&molecule,
                                                   atof(argv[2]),
                                                   grid_resolution,
                                                   grid_size,
                                                   vanderwaals,
                                                   &field);
  printf("Save Electrostatic Potential ...\n");
  SaveDXField(field, argv[5]);
  DelVoxel(&field);
  DelMolecule(&molecule);
  return 0;
}
