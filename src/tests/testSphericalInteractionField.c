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
  if(argc != 4){
    printf("\nUsage %s file.mol2 [forml charge] forcefield.txt outfield.mol2\n\n", argv[0]);
    printf("When you get the output you can load this into pymol and color by field value by using \"spectrum pc, blue_red, O_Testosterone.fiedl\"\n");
    return -1;
  }

  int formal_charge = atoi(argv[2]);

  ForceField ff;
  LoadFFParams(argv[3], &ff);

  /*PrintFFParams(ff);*/

  matrix *field;
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  printf("%s", molecule.molname);

  double steps[4][2] = {{0.00, -1.0}, {-1.0, -2.0}, {-4.0, -2.0}, {-100.0, -4.0}};

  //int probe_id[4] = {0, 32, 2, 29};
  //char ptype[4][10] = {{"N"}, {"O_Water"}, {"C"}, {"H_Water"}};

  int probe_id[2] = {0, 1};
  char ptype[4][10] = {{"O_Water"}, {"H_Water"}};

  for(int i = 0; i < 2; i++){
    initMatrix(&field);
    /*FieldCalculator(ForceField ff, MOLECULE molecule, size_t npnt, enum RADIUS_TYPE rtype, int probe_id, matrix **field);*/
    FieldCalculator(ff,
                    &molecule,
                    formal_charge,
                    500,
                    vanderwaals,
                    probe_id[i],
                    0.f,
                    field);
    matrix *spectra;
    initMatrix(&spectra);
    SphericalFieldAnalysis(molecule, field, 0.025, spectra);
    /*Normalize between 0 and 1*/
    PrintMatrix(spectra);
    DelMatrix(&spectra);

    for(int j = 0; j < 4; j++){
      printf("\t%f", AvgFieldValue(field, steps[j][0], steps[j][1]));
      AvgFieldValue(field, steps[j][0], steps[j][1]);
    }
    char *s = concat(3, ptype[i], "_", argv[4]);
    WriteMol2SphericalPoints(field, s);
    free(s);
    /*SaveCGOField(field, argv[4]);*/
    DelMatrix(&field);
  }
  printf("\n");
  /*MEPCalculator(molecule, 0.25, vanderwaals, &mep, NULL);*/
  /*SaveXPLORField(mep, argv[2]);
  DeleteMEP(&mep);*/
  DelMolecule(&molecule);
  DeleteForceField(&ff);
  return 0;
}
