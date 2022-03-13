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
    printf("\nUsage %s file.mol2 [txt param file] [grid resolution] [grid size] outfield.dx\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  printf("Load molecule %s ...\n", molecule.molname);

  int grid_resolution = atoi(argv[3]);
  int grid_size = atoi(argv[4]);

  printf("Read atom parameters...\n");
  AtomsProperty *lst;
  ReadAtomProperties(argv[2], &lst);

  printf("Calculate Generic Potential ...\n");
  VOXEL *field;
  VoxelGenericPotentialCalculator(&molecule, lst, grid_resolution, grid_size, vanderwaals, &field);
  printf("Save Generic Potential ...\n");
  SaveDXField(field, argv[5]);

  //SaveMol2VoxelField(field, argv[5]);

  DelVoxel(&field);
  DelMolecule(&molecule);
  DeleteAtomProperties(&lst);
  return 0;
}
