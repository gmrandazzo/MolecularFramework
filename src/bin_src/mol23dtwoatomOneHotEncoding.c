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
  if(argc != 7){
    printf("\nUsage %s [file.mol2] [atom_id1] [atom_id2] [sphere step] [max radius] [CSV One Hot Encoding]\n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    int atom_id1 = atoi(argv[2]);
    int atom_id2 = atoi(argv[3]);
    double sph_step = atof(argv[4]);
    double max_radius = atof(argv[5]);
    dvector *fp;
    Get3DTwoAtomSumDistanceFingerprint(molecule, atom_id1, atom_id2, sph_step, max_radius, &fp);

    /*
    dvector *fp_dist, *fp_dihedral;
    Get3DTwoAtomSumDistanceFingerprint(molecule, atom_id1, atom_id2, sph_step, max_radius, &fp_dist);

    Get3DTwoAtomBaricentreSumDihedralAngleFingerprint(molecule,
                                                      atom_id1,
                                                      atom_id2,
                                                      sph_step,
                                                      max_radius,
                                                      &fp_dihedral);
    */

    /*
    Combined version of Distance*dihedral angle...
    Get3DTwoAtomDistanceDihedralAngleCombinationFingerprint(molecule,
                                                            atom_id1,
                                                            atom_id2,
                                                            sph_step,
                                                            max_radius,
                                                            &fp);
                                                            */
    FILE *f = fopen(argv[6], "w");
    matrix *m;
    initMatrix(&m);
    FP2Matrix(fp, &m);
    PrintMatrix(m);
    for(int i = 0; i < m->row; i++){
      for(int j = 0; j < m->col-1; j++){
        fprintf(f, "%.8f,", m->data[i][j]);
      }
      fprintf(f, "%.8f\n", m->data[i][m->col-1]);
    }
    fclose(f);
    DelMatrix(&m);
    DelMolecule(&molecule);
    DelDVector(&fp);
    return 0;
  }
}
