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
#include "../fielddesc.h"

int main(int argc, char **argv)
{
  if(argc != 8){
    printf("\nUsage %s file.mol2 [formal charge] forcefield.txt [probe id in forcefield.txt] [grid resolution] [grid size] outfield.dx\n\n", argv[0]);
    return -1;
  }

  ForceField ff;
  printf("Load probe params ...\n");
  LoadFFParams(argv[3], &ff);

  /*PrintFFParams(ff);*/

  MOLECULE molecule;
  printf("Load molecule %s ...\n", molecule.molname);
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 5);

  int formal_charge = atoi(argv[2]);
  int probe_id = atoi(argv[4]);
  int grid_resolution = atoi(argv[5]);
  int grid_size = atoi(argv[6]);

  printf("Calculate MIFs ...\n");
  VOXEL *field;
  VoxelFieldCalculator(&molecule, formal_charge, grid_resolution, grid_size, vanderwaals, probe_id, ff, &field);
  printf("Save MIFs ...\n");
  SaveDXField(field, argv[7]);

  DelVoxel(&field);
  DelMolecule(&molecule);
  DeleteForceField(&ff);
  return 0;
}
