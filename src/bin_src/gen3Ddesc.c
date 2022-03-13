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
#include "../mol3Daligner.h"
#include "../mol3Dfields.h"
#include "../fielddesc.h"

int main(int argc, char **argv)
{
  if(argc != 3){
    printf("\nUsage %s file.mol2 [formal charge] forcefield.txt\n\n", argv[0]);
    return -1;
  }

  int n_esteps = 16;
  /*double esteps[8][2] = {{-10, -3.5}, {-3.5, -3.0}, {-3.0, -2.5}, {-2.5, -2.0}, {-2.0, -1.5}, {-1.5, -1.0}, {-1.0, -0.5}, {-0.5, 0}};*/
  double esteps[16][2] = {{0.00, -1.0},
                          {-2.0, -1.0},
                          {-3.0, -2.0},
                          {-4.0, -3.0},
                          {-5.0, -4.0},
                          {-6.0, -5.0},
                          {-7.0, -7.0},
                          {-8.0, -7.0},
                          {-10.0, -8.0},
                          {-14.0, -10.0},
                          {-20.0, -14.0},
                          {-30.0, -20.0},
                          {-60.0, -30.0},
                          {-100.0, -60.0},
                          {-200.0, -100.0},
                          {-600.0, -200.0}};

  ForceField ff;
  LoadFFParams(argv[3], &ff);

  matrix *Nfield;
  matrix *O_WATfield;
  matrix *C_R_field;
  matrix *C_aro_field;
  matrix *H_field;

  initMatrix(&Nfield);
  initMatrix(&O_WATfield);
  initMatrix(&C_R_field);
  initMatrix(&C_aro_field);
  initMatrix(&H_field);

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  /*random conformational rototranslation to be rotational/translational invariant*/
  RandomConformationRotation(&molecule);

  AtomAnalyzer(&molecule, 1);
  printf("%s", molecule.molname);

  int formal_charge = atoi(argv[2]);
  int n_fpnt = 200;
   /* N as NH3 Hydrogen bond donor to describe hydrogen bond acceptors in molecule */
  FieldCalculator(ff, &molecule, formal_charge, 100, vanderwaals, 34, 0.f, &Nfield);
  /*RefitFields(&Nfield, -1, 1);*/
  /* O as Hydrogen bond acceptor to describe hydrogen bond donor in molecule */
  FieldCalculator(ff, &molecule, formal_charge, 100, vanderwaals, 32, 0.f, &O_WATfield);
  /*RefitFields(&O_WATfield, -1, 1);*/
  /* C as Hydrophobic interaction with an aliphatic in molecule */
  FieldCalculator(ff, &molecule, formal_charge,  n_fpnt, vanderwaals, 36, 0.f, &C_R_field);
  /* C as Hydrophobic interaction with an aromatic in molecule  */
  FieldCalculator(ff, &molecule, formal_charge, n_fpnt, vanderwaals, 38, 0.f, &C_aro_field);
  /*RefitFields(&C_field, -1, 1);*/
  /* H a shape interaction with an hydrogen in molecule  */
  FieldCalculator(ff, &molecule, formal_charge, n_fpnt, vanderwaals, 32, 0.f, &H_field);
  /*RefitFields(&H_field, -1, 1);*/

  /* Single field descriptors */
  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f\t%f", AvgFieldValue(Nfield, esteps[j][0], esteps[j][1]),
                               DipoleFieldValue(Nfield, esteps[j][0], esteps[j][1]),
                               BaricenterFieldValue(molecule, Nfield, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f\t%f", AvgFieldValue(O_WATfield, esteps[j][0], esteps[j][1]),
                               DipoleFieldValue(O_WATfield, esteps[j][0], esteps[j][1]),
                               BaricenterFieldValue(molecule, O_WATfield, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f\t%f", AvgFieldValue(C_R_field, esteps[j][0], esteps[j][1]),
                               DipoleFieldValue(C_R_field, esteps[j][0], esteps[j][1]),
                               BaricenterFieldValue(molecule, C_R_field, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f\t%f", AvgFieldValue(C_aro_field, esteps[j][0], esteps[j][1]),
                               DipoleFieldValue(C_aro_field, esteps[j][0], esteps[j][1]),
                               BaricenterFieldValue(molecule, C_aro_field, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f\t%f", AvgFieldValue(H_field, esteps[j][0], esteps[j][1]),
                               DipoleFieldValue(H_field, esteps[j][0], esteps[j][1]),
                               BaricenterFieldValue(molecule, H_field, esteps[j][0], esteps[j][1]));
  }

  /* Double field descriptors */
  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f", DiffFieldValue(Nfield, H_field, esteps[j][0], esteps[j][1]),
                       RatioFieldValue(Nfield, H_field, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f", DiffFieldValue(O_WATfield, H_field, esteps[j][0], esteps[j][1]),
                       RatioFieldValue(O_WATfield, H_field, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f", DiffFieldValue(C_R_field, H_field, esteps[j][0], esteps[j][1]),
                       RatioFieldValue(C_R_field, H_field, esteps[j][0], esteps[j][1]));
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f\t%f", DiffFieldValue(C_aro_field, H_field, esteps[j][0], esteps[j][1]),
                       RatioFieldValue(C_aro_field, H_field, esteps[j][0], esteps[j][1]));
  }

  printf("\n");

  DelMolecule(&molecule);
  DelMatrix(&Nfield);
  DelMatrix(&O_WATfield);
  DelMatrix(&C_R_field);
  DelMatrix(&C_aro_field);
  DelMatrix(&H_field);
  DeleteForceField(&ff);
  return 0;
}
