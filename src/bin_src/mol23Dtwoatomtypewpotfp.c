/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <scientific.h>
#include "../geomdesc.h"
#include "../molecule.h"

int main(int argc, char **argv)
{
  if(argc != 7){
    printf("\nUsage %s [file.mol2] [txt param file] [atom_id1] [atom_id2] [sphere step] [max radius]\n\n", argv[0]);
    return -1;
  }
  else{

    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    int atom_id1 = atoi(argv[3]);
    int atom_id2 = atoi(argv[4]);
    double sph_step = atof(argv[5]);
    double max_radius = atof(argv[6]);
    int i;
    dvector *fp_dist, *fp_alpha_a1, *fp_alpha_a2, *fp_dihedral_a1, *fp_dihedral_a2;

    AtomsProperty *a_weights;
    ReadAtomProperties(argv[2], &a_weights);

    Get3DTwoAtomSumPotentialFingerprint(molecule, atom_id1, atom_id2, sph_step, max_radius, a_weights, &fp_dist);
    //Get3DTwoAtomSumDistanceFingerprint(molecule, atom_id1, atom_id2, sph_step, max_radius, &fp_dist);

    GetAngleSumFingerprint(molecule, atom_id1, &fp_alpha_a1);
    GetTorsionSumFingerprint(molecule, atom_id1, &fp_dihedral_a1);

    GetAngleSumFingerprint(molecule, atom_id2, &fp_alpha_a2);
    GetTorsionSumFingerprint(molecule, atom_id2, &fp_dihedral_a2);

    printf("%s_%d%s->%d%s,", molecule.molname,
                             atom_id1,
                             molecule.atoms[atom_id1].type,
                             atom_id2,
                             molecule.atoms[atom_id2].type);

    for(i = 0; i < fp_dist->size; i++){
      printf("%.8f,", fp_dist->data[i]);
    }

    for(i = 0; i < fp_alpha_a1->size; i++){
      printf("%.8f,", fp_alpha_a1->data[i]+fp_alpha_a2->data[i]);
    }

    for(i = 0; i < fp_dihedral_a1->size-1; i++){
      printf("%.8f,", fp_dihedral_a1->data[i]+fp_dihedral_a2->data[i]);
    }

    i = fp_dihedral_a1->size-1;
    printf("%.8f\n", fp_dihedral_a1->data[i]+fp_dihedral_a2->data[i]);

    DelMolecule(&molecule);
    DelDVector(&fp_dist);
    DelDVector(&fp_alpha_a1);
    DelDVector(&fp_alpha_a2);
    DelDVector(&fp_dihedral_a1);
    DelDVector(&fp_dihedral_a2);
    DeleteAtomProperties(&a_weights);
    return 0;
  }
}
