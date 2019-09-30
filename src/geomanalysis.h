/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef GEOMANALYSIS_H
#define GEOMANALYSIS_H
#include <iostream>
#include <string>
#include <vector>
#include "molecule.h"

extern "C" {
    #include <scientific.h>
}


using namespace std;

class GeomAnalysis
{
public:
  GeomAnalysis(){}
  double SimplexAlign3DConformations(MOLECULE *from, MOLECULE *to);

private:
  matrix *_Aconf_;
  matrix *_Bconf_;
  /* Utilised by NelderMead to minimize something... */
  double func(dvector *x);
  void ComputeCoordinates_(matrix *A_T, matrix *R, dvector *t, matrix **A_aligned_T);
  void ConvertToRandt_v2_(dvector *x, matrix *R, dvector *t);
  void ConvertToRandt_v1_(dvector *x, matrix *R, dvector *t);

};

#endif
