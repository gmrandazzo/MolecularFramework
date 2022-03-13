/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <scientific.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../mol3Daligner.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf(" Align a 3D molecule A to a molecule B\n");
    printf("\nUsage %s align_A.mol2 in_B.mol2 [map file A->B]\n\n", argv[0]);
    return -1;
  }

  MOLECULE A, B;
  NewMOL2Molecule(&A, argv[1]);
  AtomAnalyzer(&A, 1);
  NewMOL2Molecule(&B, argv[2]);
  AtomAnalyzer(&B, 1);

  uivector *aidA, *aidB;
  initUIVector(&aidA);
  initUIVector(&aidB);

  int j;
  char l[2048], *tl;
  FILE *fi;
  if((fi = fopen (argv[3], "r")) == NULL){
     printf ("File %s not found!!!", argv[3]);
     return 0;
  }

  fseek (fi, 0, SEEK_SET);
  while(fgets(l, LINE_MAX, fi) != NULL){
    // Tokenize
    for(tl = strtok(Trim(l), " \t"), j = 0; tl; tl = strtok (NULL, " \t"), j++){
      if(j == 0){
        UIVectorAppend(&aidA, atoi(tl));
      }
      else if(j == 1){
        UIVectorAppend(&aidB, atoi(tl));
      }
      else{
        continue;
      }
    }
  }
  fclose(fi);


  double rmsd = Align3DPharmacophore(A, B, aidA, aidB);
  printf("RMSD: %.3f\n", rmsd);

  SaveMol2Molecule(A,  argv[1]);
  DelUIVector(&aidA);
  DelUIVector(&aidB);
  DelMolecule(&A);
  DelMolecule(&B);
  return 0;
}
