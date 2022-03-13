/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <math.h>
#include "../periodic_table.h"
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Dfields.h"
#include "../fielddesc.h"


#include <stdio.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Daligner.h"
#include "../mol3Dfields.h"
#include "../fielddesc.h"

int main(int argc, char **argv)
{
  if(argc != 5){
    printf("\nUsage %s molecule_1.mol2 molecule_2.mol2 [grid resolution] [grid size]\n\n", argv[0]);
    return -1;
  }


  MOLECULE molecule1;
  NewMOL2Molecule(&molecule1, argv[1]);
  AtomAnalyzer(&molecule1, 1);
  
  MOLECULE molecule2;
  NewMOL2Molecule(&molecule2, argv[2]);
  AtomAnalyzer(&molecule2, 1);

  double rmsd = Align3DShapes(molecule1, molecule2);
  //double rmsd = Align3DOnVDWShapes(molecule1, molecule2);
  
  int grid_resolution = atoi(argv[3]);
  int grid_size = atoi(argv[4]);

  VOXEL *v1;
  VOXEL *v2;
  VoxelElectrostaticPotentialCalculator(&molecule1, grid_resolution, grid_size, vanderwaals, &v1);
  VoxelElectrostaticPotentialCalculator(&molecule2, grid_resolution, grid_size, vanderwaals, &v2);

  double norm = 0.f;
  for(int i = 0; i < v1->nx; i++){
    for(int j = 0; j < v1->ny; j++){
      for(int k = 0; k < v1->nz; k++){
        norm +=  square(v1->pnt[i][j][k] - v2->pnt[i][j][k]);
      }
    }
  }

  norm = sqrt(norm);
  printf("%s,%f,%.3f\n", molecule2.molname, norm, rmsd);

  DelVoxel(&v1);
  DelVoxel(&v2);
  DelMolecule(&molecule1);
  DelMolecule(&molecule2);
  return 0;
}
