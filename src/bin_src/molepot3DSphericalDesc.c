/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <scientific.h>
#include "../periodic_table.h"
#include "../molecule.h"
#include "../mol3Dfields.h"
#include "../mol3Daligner.h"
#include "../fielddesc.h"

int main(int argc, char **argv)
{
  if(argc != 6){
    printf("\nUsage %s file.mol2 [n spherical points] [scale min] [scale max] [n rotations]\n\n", argv[0]);
    return -1;
  }

  int npnt = atoi(argv[2]);
  double x, min = atof(argv[3]), max = atof(argv[4]);
  int nrot = atoi(argv[5]);

  int n_esteps = 10;
  double dx = (max-min)/(double)n_esteps;

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);

  dvector *freqfv, *avgfv, *dipfv, *barfv;
  NewDVector(&freqfv, n_esteps);
  NewDVector(&avgfv, n_esteps);
  NewDVector(&dipfv, n_esteps);
  NewDVector(&barfv, n_esteps);

  for(int i = 0; i < nrot; i++){
    /*random conformational rototranslation to be rotational/translational invariant*/
    RandomConformationRotation(&molecule);
    matrix *field;
    initMatrix(&field);
    SphericalElectrostaticPotentialCalculator(&molecule, npnt, vanderwaals, field);
    /*Voxel2Matrix(v, &field);*/
    //Calculate Histograms
    x = min;
    for(int j = 0; j < n_esteps; j++){
      freqfv->data[j] += (double)FreqFieldValue(field, x, x+dx);
      avgfv->data[j] += (double)AvgFieldValue(field, x, x+dx);
      dipfv->data[j] += (double)DipoleFieldValue(field, x, x+dx);
      barfv->data[j] += (double)BaricenterFieldValue(molecule, field, x, x+dx);
      x+=dx;
    }
    //DelVoxel(&v);
    DelMatrix(&field);
  }

  for(int j = 0; j < n_esteps; j++){
    freqfv->data[j] /= (double)nrot;
    avgfv->data[j] /= (double)nrot;
    dipfv->data[j] /= (double)nrot;
    barfv->data[j] /= (double)nrot;
  }

  // Output
  printf("%s", molecule.molname);
  //Histogram
  for(int j = 0; j < n_esteps; j++){
    printf("\t%f", freqfv->data[j]);
  }

  /*
  for(int j = 0; j < n_esteps; j++){
    printf("\t%f", avgfv->data[j]);
  }

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f", dipfv->data[j]);
  }
  */

  for(int j = 0; j < n_esteps; j++){
    printf("\t%f", barfv->data[j]);
  }
  printf("\n");

  DelDVector(&freqfv);
  DelDVector(&avgfv);
  DelDVector(&dipfv);
  DelDVector(&barfv);
  DelMolecule(&molecule);
  return 0;
}
