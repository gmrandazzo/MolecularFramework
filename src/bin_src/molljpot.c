/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../periodic_table.h"
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Dfields.h"

int main(int argc, char **argv)
{
  if(argc != 6){
    printf("\nUsage %s file.mol2 [epsilon txt param file] [sigma txt param file] [grid resolution] [grid size] outfield.dx\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  printf("Load molecule %s ...\n", molecule.molname);
  NewMOL2Molecule(&molecule, argv[1]);

  printf("Read epsilon atom parameters...\n");
  AtomsProperty *epsilon;
  ReadAtomProperties(argv[2], &epsilon);

  AtomsProperty *sigma;
  ReadAtomProperties(argv[3], &sigma);

  int grid_resolution = atoi(argv[4]);
  int grid_size = atoi(argv[5]);

  printf("Calculate LJ Potential ...\n");
  VOXEL *field;
  VoxelLJPotentialCalculator(&molecule, epsilon, sigma, grid_resolution, grid_size, vanderwaals, &field);

  printf("Save Generic Potential ...\n");
  SaveDXField(field, argv[6]);

  DelVoxel(&field);
  DelMolecule(&molecule);
  DeleteAtomProperties(&epsilon);
  DeleteAtomProperties(&sigma);
  return 0;
}
