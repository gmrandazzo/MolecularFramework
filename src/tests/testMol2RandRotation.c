/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <math.h>
#include "../molecule.h"
#include "../mol3Daligner.h"
#include <scientific.h>

int main(int argc, char **argv)
{
  if(argc == 3){
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    RandomConformationRotation(&molecule);
    /*
    matrix *A;
    matrix *R;
    dvector *t;
    uivector *ids;
    NewMatrix(&R, 3, 3);
    NewDVector(&t, 3);
    int nsmpl = (int)ceil(molecule.n_atoms*0.1);
    NewUIVector(&ids, nsmpl);
    printf("%d\n", nsmpl);
    NewMatrix(&A, nsmpl, 3);

    for(int i = 0; i < nsmpl; i++){
      ids->data[i] = randInt(0, molecule.n_atoms);
      printf("%zu/%d\n", ids->data[i], molecule.n_atoms);
      A->data[i][0] = molecule.atoms[ids->data[i]].coord.x;
      A->data[i][1] = molecule.atoms[ids->data[i]].coord.y;
      A->data[i][2] = molecule.atoms[ids->data[i]].coord.z;
    }

    Random3DMatrixRotationEulerMethod(A, R, t);
    RotoTranslateMolecule(molecule, R, t);
    puts("Rotation result sample-molecule");
    for(int i = 0; i < nsmpl; i++){
      printf("x: %f y: %f z: %f\n", A->data[i][0]-molecule.atoms[ids->data[i]].coord.x,
                                    A->data[i][1]-molecule.atoms[ids->data[i]].coord.y,
                                    A->data[i][2]-molecule.atoms[ids->data[i]].coord.z);
    }



    DelMatrix(&A);
    DelUIVector(&ids);
    DelMatrix(&R);
    DelDVector(&t);
    */
    SaveMol2Molecule(molecule, argv[2]);
    DelMolecule(&molecule);
    return 0;
  }
  else{
    printf("\nUsage %s infile.mol2 outfile.mol2 \n\n", argv[0]);
    return -1;
  }
}
