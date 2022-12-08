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
  if(argc != 4){
    printf("\nUsage %s file.mol2 [N voxel points (INT)] [Angstrom grid size (INT)]\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  int grid_resolution = atoi(argv[2]);
  int grid_size = atoi(argv[3]);

  VOXEL *field;
  VoxelElectrostaticPotentialCalculator(&molecule, grid_resolution, grid_size, vanderwaals, &field);
    
  printf("%s", molecule.molname);
  for(int k = 0; k < field->nz; k++){
    for(int i = 0; i < field->nx; i++){
      for(int j = 0; j < field->ny; j++){
        printf(",%f", field->pnt[i][j][k]);
      }
    }
  }
  printf("\n");
  DelVoxel(&field);
  DelMolecule(&molecule);
  return 0;
}
