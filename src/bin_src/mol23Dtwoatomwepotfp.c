/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../molecule.h"
#include "../geomdesc.h"


int main(int argc, char **argv)
{
  if(argc != 6){
    printf("\nUsage %s [file.mol2] [atom_id1] [atom_id2] [sphere step] [max radius]\n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    int atom_id1 = atoi(argv[2]);
    int atom_id2 = atoi(argv[3]);
    double sph_step = atof(argv[4]);
    double max_radius = atof(argv[5]);

    dvector *fp_epot, *fp_wdepot;

    Get3DEPotDihedralAngleWeightedTwoAtomFingerprint(molecule,
                                                     atom_id1,
                                                     atom_id2,
                                                     sph_step,
                                                     max_radius,
                                                     &fp_wdepot);

    Get3DTwoAtomEpotFingerprint(molecule,
                                atom_id1,
                                atom_id2,
                                sph_step,
                                max_radius,
                                &fp_epot);

    double dst = GetDistance_(molecule.atoms[atom_id1].coord, molecule.atoms[atom_id2].coord);
    double angle, dihedral;
    Get3DAngle(molecule, atom_id1, atom_id2, &angle, NULL);
    Get3DDihedralAngle(molecule, atom_id1, atom_id2, &dihedral, NULL);

    if(FLOAT_EQ(angle, -9999.0, 1e-3)){
      /*No angle possible*/
      angle = -9999.0;
    }
    else{
      angle = fabs(Rad2Grad(angle));
    }

    if(FLOAT_EQ(dihedral, -9999.0, 1e-3)){
      /* No dihedral possible */
      dihedral = -9999.0;
    }
    else{
      dihedral = fabs(Rad2Grad(dihedral));
    }

    printf("%s_%d%s->%d%s,%f,%f,%f,", molecule.molname,
                                      atom_id1,
                                      molecule.atoms[atom_id1].type,
                                      atom_id2,
                                      molecule.atoms[atom_id2].type,
                                      dst,
                                      angle,
                                      dihedral);
    /*
    for(int i = 0; i < fp_epot->size-1; i++){
      printf("%.8f,%.8f,%.8f,", fp_epot->data[i], fp_waepot->data[i], fp_wdepot->data[i]);
    }
    printf("%.8f,%.8f,%.8f\n", fp_epot->data[fp_epot->size-1],
                               fp_waepot->data[fp_epot->size-1],
                               fp_wdepot->data[fp_epot->size-1]);
    */

    for(int i = 0; i < fp_epot->size-1; i++){
      printf("%.8f,%.8f,", fp_epot->data[i], fp_wdepot->data[i]);
    }

    printf("%.8f,%.8f\n", fp_epot->data[fp_epot->size-1],
                          fp_wdepot->data[fp_epot->size-1]);
    DelMolecule(&molecule);
    DelDVector(&fp_wdepot);
    /*DelDVector(&fp_waepot);*/
    DelDVector(&fp_epot);
    return 0;
  }
}
