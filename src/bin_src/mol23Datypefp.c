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
#include "../atomanalysis.h"


int main(int argc, char **argv)
{
  if(argc <= 3 || argc >= 6){
    printf("\nUsage for a specific atom: %s [file.mol2] [atom_id] [sphere step] [max radius]\n\n", argv[0]);
    printf("\nUsage for all atoms      : %s [file.mol2] [sphere step] [max radius]\n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    AtomAnalyzer(&molecule, 3); /* Necessary for hyprid type */
    int atom_id = 0;
    double sph_step = 0.f;
    double max_radius = 0.f;
    dvector *fp, *fp_alpha, *fp_dihedral;/*, *fp_bari;*/
    if(argc == 5){
      atom_id = atoi(argv[2]);
      sph_step = atof(argv[3]);
      max_radius = atof(argv[4]);
      Get3DAtomSumDistanceFingerprint(molecule, atom_id, sph_step, max_radius, &fp);
      GetAngleSumFingerprint(molecule, atom_id, &fp_alpha);
      GetTorsionSumFingerprint(molecule, atom_id, &fp_dihedral);
      /*Get3DBaricentreSumAngleAtomFingerprint(molecule, atom_id, sph_step, max_radius, &fp_bari);*/
      printf("%s_%d_%s_", molecule.molname, atom_id+1, molecule.atoms[atom_id].asymbl);
      int t = molecule.atoms[atom_id].ainfo.hybrid;
      switch(t){
        case 1:
            printf("S,");
            break;
        case 2:
            printf("SP,");
            break;
        case 3:
            printf("SP2,");
            break;
        case 4:
            printf("SP3,");
            break;
        case 5:
            printf("SP3D,");
            break;
        case 6:
            printf("SP3D2,");
            break;
        case 7:
            printf("AR,");
            break;
        case 8:
            printf("OTHER,");
            break;
        default:
            printf("UNSPECIFIED,");
      }
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
      /*
      for(int i = 0; i < fp_bari->size-1; i++){
        printf("%.8f,", fp_bari->data[i]);
      }
      printf("%.8f\n", fp_bari->data[fp_bari->size-1]);
      */
      DelDVector(&fp);
      DelDVector(&fp_alpha);
      DelDVector(&fp_dihedral);
      /*DelDVector(&fp_bari);*/
    }
    else{
      sph_step = atof(argv[2]);
      max_radius = atof(argv[3]);
      for(atom_id = 0; atom_id < molecule.n_atoms; atom_id++){
        Get3DAtomSumDistanceFingerprint(molecule, atom_id, sph_step, max_radius, &fp);
        GetAngleSumFingerprint(molecule, atom_id, &fp_alpha);
        GetTorsionSumFingerprint(molecule, atom_id, &fp_dihedral);
        /*Get3DBaricentreSumAngleAtomFingerprint(molecule, atom_id, sph_step, max_radius, &fp_bari);*/
        printf("%s_%d_%s_", molecule.molname, atom_id+1, molecule.atoms[atom_id].asymbl);
        int t = molecule.atoms[atom_id].ainfo.hybrid;
        switch(t){
            case 1:
                printf("S,");
                break;
            case 2:
                printf("SP,");
                break;
            case 3:
                printf("SP2,");
                break;
            case 4:
                printf("SP3,");
                break;
            case 5:
                printf("SP3D,");
                break;
            case 6:
                printf("SP3D2,");
                break;
            case 7:
                printf("AR,");
                break;
            case 8:
                printf("OTHER,");
                break;
            default:
                printf("UNSPECIFIED,");
        }
        
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
        /*
        for(int i = 0; i < fp_bari->size-1; i++){
          printf("%.8f,", fp_bari->data[i]);
        }
        printf("%.8f\n", fp_bari->data[fp_bari->size-1]);
        */
        DelDVector(&fp);
        DelDVector(&fp_alpha);
        DelDVector(&fp_dihedral);
        /*DelDVector(&fp_bari);*/
      }
    }
    DelMolecule(&molecule);
  }
  return 0;
}
