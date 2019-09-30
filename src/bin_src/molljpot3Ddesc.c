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
  if(argc != 9){
    printf("\nUsage %s file.mol2 [epsilon txt param file] [sigma txt param file] [grid resolution] [grid size] [scale min] [scale max] [n rot]\n\n", argv[0]);
    return -1;
  }

  int nrot = atoi(argv[8]);
  int n_esteps = 10;
  double x, min = atof(argv[5]), max = atof(argv[6]);
  double dx = (max-min)/(double)n_esteps;

  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);

  /*random conformational rototranslation to be rotational/translational invariant*/
  RandomConformationRotation(&molecule);

  AtomsProperty *epsilon;
  ReadAtomProperties(argv[2], &epsilon);

  AtomsProperty *sigma;
  ReadAtomProperties(argv[3], &sigma);

  int grid_resolution = atoi(argv[4]);
  int grid_size = atoi(argv[5]);

  dvector *freqfv, *avgfv, *dipfv, *barfv;
  NewDVector(&freqfv, n_esteps);
  NewDVector(&avgfv, n_esteps);
  NewDVector(&dipfv, n_esteps);
  NewDVector(&barfv, n_esteps);

  for(int i = 0; i < nrot; i++){
    /*random conformational rototranslation to be rotational/translational invariant*/
    RandomConformationRotation(&molecule);
    VOXEL *v;
    VoxelLJPotentialCalculator(&molecule, epsilon, sigma, grid_resolution, grid_size, vanderwaals, &v);
    matrix *field;
    initMatrix(&field);
    Voxel2Matrix(v, &field);
    //Calculate Histograms
    x = min;
    for(int j = 0; j < n_esteps; j++){
      freqfv->data[j] += (double)FreqFieldValue(field, x, x+dx);
      avgfv->data[j] += (double)AvgFieldValue(field, x, x+dx);
      dipfv->data[j] += (double)DipoleFieldValue(field, x, x+dx);
      barfv->data[j] += (double)BaricenterFieldValue(molecule, field, x, x+dx);
      x+=dx;
    }

    /*
    matrix *m;
    initMatrix(&m);
    Voxel2Matrix(v, &m);
    matrix *spectra;
    initMatrix(&spectra);
    SphericalFieldAnalysis(molecule, m, 0.025, &spectra);
    PrintMatrix(spectra);
    DelMatrix(&m);
    DelMatrix(&spectra);
    SaveDXField(v, "test.dx");
    SaveMol2VoxelField(v,"test.mol2");
    */

    DelVoxel(&v);
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
  DeleteAtomProperties(&sigma);
  DeleteAtomProperties(&epsilon);
  return 0;
}
