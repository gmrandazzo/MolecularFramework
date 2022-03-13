/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../wfingerprint.h"


int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage\n");
    printf("  >> Init Table: %s 0 <weightable.txt> <file1.mol2> <file2.mol2> <file3.mol2> <...> <fileN.mol2>\n", argv[0]);
    printf("  >> Calculate Sum Property: %s 1 <weightable > <file1.mol2>\n", argv[0]);
    printf("  >> Calculate Product Property: %s 2 <weightable> <file1.mol2>\n", argv[0]);
    printf("  >> Write Atom Contribution: %s 3 <weightable> <file1.mol2> <file1_output_contribution.txt>\n", argv[0]);
    printf("  >> Get Atom Fragment Fingerprint: %s 4 <weightable> <file1.mol2>\n\n", argv[0]);
    printf("%s was written by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>\n\n", argv[0]);
    return -1;
  }
  else{
    int what = atoi(argv[1]);
    int i;
    double property;
    MOLECULE molecule;
    if(what == 0){
      for(i = 3; i < argc; i++){
        NewSDFMolecule(&molecule, argv[i]);
        AtomAnalyzer(&molecule, 1);
        //NewMOL2Molecule(&molecule, argv[i]);
        if(molecule.n_atoms > 0 && molecule.n_bonds > 0)
          InitWeightTable(argv[2], molecule);
        DelMolecule(&molecule);
      }
    }
    else if(what == 1){
      NewSDFMolecule(&molecule, argv[3]);
      AtomAnalyzer(&molecule, 1);
      //NewMOL2Molecule(&molecule, argv[3]);
      CalcPropertFromSumAtomFrag(molecule, argv[2], &property);
      printf("%s,%f\n", molecule.molname, property);
      DelMolecule(&molecule);
    }
    else if(what == 2){
      NewSDFMolecule(&molecule, argv[3]);
      AtomAnalyzer(&molecule, 1);
      //NewMOL2Molecule(&molecule, argv[3]);
      CalcPropertFromProductAtomFrag(molecule, argv[2], &property);
      printf("%s,%f\n", molecule.molname, property);
      DelMolecule(&molecule);
    }
    else if(what == 3){
      if(argc == 5){
        NewSDFMolecule(&molecule, argv[3]);
        AtomAnalyzer(&molecule, 1);
        //NewMOL2Molecule(&molecule, argv[3]);
        WriteAtomFragWeight(molecule, argv[2], argv[4]);
        DelMolecule(&molecule);
      }
      else{
        printf("Write Atom Contribution: %s 3 <weightable> <file1.mol2> <file1_output_contribution.txt>\n\n", argv[0]);
      }
    }
    else if(what == 4){
      if(argc == 4){
        double *fp = NULL;
        int fpsize;
        NewSDFMolecule(&molecule, argv[3]);
        AtomAnalyzer(&molecule, 1);
        //NewMOL2Molecule(&molecule, argv[3]);
        GetSumAtomFragFP(molecule, argv[2], &fp, &fpsize);
        printf("%s,", molecule.molname);
        for(int i = 0; i < fpsize-1; i++)
          printf("%f,", fp[i]);
        printf("%f\n", fp[fpsize-1]);
        free(fp);
        DelMolecule(&molecule);
      }
      else{
        printf("Write AtomFragment Fingerprint: %s 4 <weightable> <file1.mol2>\n\n", argv[0]);
      }
    }
    else{
      printf("Unknown option selected!\n");
    }
    return 0;
  }
}
