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
  if(argc != 5){
    printf("\nUsage %s file.mol2 [N voxel points (INT)] [Angstrom grid size (INT)] outfield.dx\n\n", argv[0]);
    return -1;
  }

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  printf("Load molecule %s ...\n", molecule.molname);

  int grid_resolution = atoi(argv[2]);
  int grid_size = atoi(argv[3]);

  printf("Calculate Electrostatic Potential ...\n");
  VOXEL *field;
  VoxelElectrostaticPotentialCalculator(&molecule, grid_resolution, grid_size, vanderwaals, &field);
  printf("Save Electrostatic Potential ...\n");
  SaveDXField(field, argv[4]);
  /* for CPCA
  for(int k = 0; k < field->nz; k++){
    printf("Block %d\n", k);
    for(int i = 0; i < field->nx; i++){
      for(int j = 0; j < field->ny-1; j++){
        printf("%f,", field->pnt[i][j][k]);
      }
      printf("%f\n", field->pnt[i][field->ny-1][k]);
    }
  }*/
  DelVoxel(&field);
  DelMolecule(&molecule);
  return 0;
}
