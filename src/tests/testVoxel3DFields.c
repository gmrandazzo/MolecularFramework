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
  if(argc != 5){
    printf("\nUsage %s file.mol2 [formal charge] forcefield.txt [n points] outfield.dx\n\n", argv[0]);
    printf("When you get the DX output you can load this into VMD and color by isovalue and so on.\n");
    return -1;
  }

  int formal_charge = atoi(argv[2]);

  ForceField ff;
  LoadFFParams(argv[3], &ff);

  /*PrintFFParams(ff);*/

  VOXEL *field;
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);

  printf("%s", molecule.molname);

  int probe_id[2] = {0, 1};
  char ptype[4][10] = {{"O_Water"}, {"H_Water"}};

  int npnts = atoi(argv[4]);

  for(int i = 0; i < 2; i++){
    VoxelFieldCalculator(&molecule, formal_charge, npnts, 8, vanderwaals, probe_id[i], ff, &field);
    char *s = concat(3, ptype[i], "_", argv[5]);
    SaveDXField(field, s);
    //SaveCubeFile(molecule, field, s);
    free(s);
    DelVoxel(&field);
  }

  DelMolecule(&molecule);
  DeleteForceField(&ff);
  return 0;
}
