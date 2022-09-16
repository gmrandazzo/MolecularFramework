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
  if(argc != 7){
    printf("\nUsage %s protein.mol2 ligand.mol2 [N voxel points (INT)] [Angstrom grid size (INT)] protein.dx ligand.dx\n\n", argv[0]);
    return -1;
  }

  MOLECULE protein;
  MOLECULE ligand;
  NewMOL2Molecule(&protein, argv[1]);
  printf("Load Protein %s ...\n", protein.molname);
  NewMOL2Molecule(&ligand, argv[2]);
  printf("Load Ligand %s ...\n", ligand.molname);

  int grid_resolution = atoi(argv[3]);
  int grid_size = atoi(argv[4]);

  printf("Calculate Electrostatic Potential ...\n");
  VOXEL *ligand_field;
  VOXEL *pocket_field;
  ProteinLigandVoxelElectrostaticPotentialCalculator(&protein,
                                                     &ligand,
                                                     grid_resolution,
                                                     grid_size,
                                                     vanderwaals,
                                                     &pocket_field,
                                                     &ligand_field);
  printf("Save Electrostatic Potential ...\n");
  SaveDXField(pocket_field, argv[5]);
  SaveDXField(ligand_field, argv[6]);

  DelVoxel(&pocket_field);
  DelVoxel(&ligand_field);
  DelMolecule(&ligand);
  DelMolecule(&protein);
  return 0;
}
