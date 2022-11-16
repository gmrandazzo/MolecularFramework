/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <math.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Dfields.h"

int main(int argc, char **argv)
{
  if(argc != 5){
    printf("\nUsage %s file1.mol2 file2.mol2 [N voxel points (INT)] [Angstrom grid size (INT)]\n\n", argv[0]);
    return -1;
  }

  MOLECULE m1;
  NewMOL2Molecule(&m1, argv[1]);

  MOLECULE m2;
  NewMOL2Molecule(&m2, argv[2]);
  //printf("Molecule %s and %s loaded\n", m1.molname, m2.molname);

  int grid_resolution = atoi(argv[3]);
  int grid_size = atoi(argv[4]);

  //printf("Calculate Electrostatic Potential ...\n");
  VOXEL *f1;
  VoxelElectrostaticPotentialCalculator(&m1, grid_resolution, grid_size, vanderwaals, &f1);

  VOXEL *f2;
  VoxelElectrostaticPotentialCalculator(&m2, grid_resolution, grid_size, vanderwaals, &f2);

  double diff = 0.f;
  for(int i = 0; i < f1->nx; i++){
    for(int j = 0; j < f1->ny; j++){
      for(int k = 0; k < f1->nz; k++){
        diff += pow(f1->pnt[i][j][k]-f2->pnt[i][j][k], 2);
      }
    }
  }

  printf("%s-%s: %f\n", m1.molname, m2.molname, diff);
  DelVoxel(&f2);
  DelVoxel(&f1);
  DelMolecule(&m1);
  DelMolecule(&m2);
  return 0;
}
