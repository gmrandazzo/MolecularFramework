/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include "../molecule.h"
#include "../geomdesc.h"


int main(int argc, char **argv)
{
  if(argc != 5){
    printf("\nUsage %s [file.mol2] [atom_id] [sphere step] [max radius]\n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    int atom_id = atoi(argv[2]);
    double sph_step = atof(argv[3]);
    double max_radius = atof(argv[4]);
    dvector *fp, *fp_alpha, *fp_dihedral;
    Get3DAtomSumDistanceFingerprint(molecule, atom_id, sph_step, max_radius, &fp);
    GetAngleSumFingerprint(molecule, atom_id, &fp_alpha);
    GetTorsionSumFingerprint(molecule, atom_id, &fp_dihedral);
    //Get3DBaricentreSumAngleAtomFingerprint(molecule, atom_id, sph_step, max_radius, &fp_alpha);
    printf("%s_%d,", molecule.molname, atom_id);
    for(int i = 0; i < fp->size; i++){
      printf("%.8f,", fp->data[i]);
    }
    for(int i = 0; i < fp_alpha->size; i++){
      printf("%.8f,", fp_alpha->data[i]);
    }
    for(int i = 0; i < fp_dihedral->size-1; i++){
      printf("%.8f,", fp_dihedral->data[i]);
    }
    printf("%.8f\n", fp_dihedral->data[fp_dihedral->size-1]);

    DelMolecule(&molecule);
    DelDVector(&fp);
    DelDVector(&fp_alpha);
    DelDVector(&fp_dihedral);
    return 0;
  }
}
