/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1 /*afprintf problem*/
#endif

#include <stdio.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Daligner.h"
#include "../mol3Dfields.h"

int main(int argc, char **argv)
{
  if(argc == 6){
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);

    int nrot = atoi(argv[2]);
    int grid_size = atoi(argv[3]); /*cube dimension*/
    int grid_resolution = atoi(argv[4]); /*number voxel point*/
    if(nrot > 0){
      for(int i = 0; i < nrot; i++){
        RandomConformationRotation(&molecule);
        /*SaveMol2Molecule(molecule, "out.mol2");*/
        VOXEL *field;
        VoxelElectrostaticPotentialCalculator(&molecule, grid_resolution, grid_size, vanderwaals, &field);
        char* new_name;
        asprintf(&new_name, "%s/%s.%d.epot.dx", argv[5], molecule.molname, i);
        SaveDXField(field, new_name);
        DelVoxel(&field);
        free(new_name);
      }
    }
    else{
      VOXEL *field;
      VoxelElectrostaticPotentialCalculator(&molecule, grid_resolution, grid_size, vanderwaals, &field);
      char* new_name;
      asprintf(&new_name, "%s/%s.%d.epot.dx", argv[5], molecule.molname, 0);
      SaveDXField(field, new_name);
      DelVoxel(&field);
      free(new_name);
    }
    DelMolecule(&molecule);
    return 0;
  }
  else{
    printf("\nUsage %s [mol2 input] [n. rotations] [grid size] [grid resolution] [output directory] \n\n", argv[0]);
    return -1;
  }
}
