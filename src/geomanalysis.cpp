/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "geomanalysis.h"
#include "mol3Daligner.h"
#include <scientific.h>
#include <cmath>

void GeomAnalysis::ComputeCoordinates_(matrix *A_T, matrix *R, dvector *t, matrix **A_aligned_T)
{
  ResizeMatrix(A_aligned_T, A_T->row, A_T->col);

  MatrixDotProduct(R, A_T, (*A_aligned_T));

  for(int i = 0; i < A_T->col; i++){
    for(int j = 0; j < A_T->row; j++)
      (*A_aligned_T)->data[j][i] += t->data[j];
  }
}


void GeomAnalysis::ConvertToRandt_v1_(dvector *x, matrix *R, dvector *t)
{
  size_t i, j, k;
  k = 0;
  for(i = 0; i < R->row; i++){
    for(j = 0; j < R->col; j++){
      R->data[i][j] = x->data[k];
      k++;
    }
  }

  for(i = 0; i < t->size; i++){
    t->data[i] = x->data[k];
    k++;
  }
}

void GeomAnalysis::ConvertToRandt_v2_(dvector *x, matrix *R, dvector *t)
{
  /*x size = 6 (3 angles + 3 translation values) */

  double A, B, C;
  A = x->data[0];
  B = x->data[1];
  C = x->data[2];

  if(A < 0.f || A > 90.f) A = 90.f;
  if(B < 0.f || B > 90.f) B = 90.f;
  if(C < 0.f || C > 90.f) C = 90.f;

  printf("A: %f B: %f C: %f\n", A, B, C);

  R->data[0][0] = cos(C)*cos(B);  R->data[0][1] = -1*sin(C); R->data[2][2] = cos(C)*sin(B);
  R->data[1][0] = cos(A)*sin(C)*cos(B)+sin(A)*sin(B); R->data[1][1] = cos(A)*cos(C); R->data[2][2] = cos(A)*sin(C)*sin(B)-sin(A)*cos(B);
  R->data[2][0] = sin(A)*sin(C)*cos(B)-cos(A)*sin(B); R->data[2][1] = sin(A)*cos(C); R->data[2][2] = sin(A)*sin(C)*sin(B)+cos(A)*cos(B);

  t->data[0] = x->data[3];
  t->data[1] = x->data[4];
  t->data[2] = x->data[5];
}

double GeomAnalysis::func(dvector *x)
{
  size_t i, j;
  matrix *R_;
  dvector *t_;
  NewMatrix(&R_, 3, 3);
  NewDVector(&t_, 3);

  PrintDVector(x);
  //ConvertToRandt_v1_(x, R_, t_);
  ConvertToRandt_v2_(x, R_, t_);

  /* Now change the coordinates to m1 with R and t
   * in order to make effective the alignment to m2 */
  matrix *A_T;
  NewMatrix(&A_T, _Aconf_->col, _Aconf_->row);
  MatrixTranspose(_Aconf_, A_T);

  matrix *A_aligned_T_tmp;
  initMatrix(&A_aligned_T_tmp);

  ComputeCoordinates_(A_T, R_, t_, &A_aligned_T_tmp);

  double rmsd = RMSD(A_aligned_T_tmp, _Bconf_);

  printf("%f\n", rmsd);

  matrix *c1, *c2;
  initMatrix(&c1);
  initMatrix(&c2);
  MatrixCovariance(A_T, &c1);
  MatrixCovariance(A_aligned_T_tmp, &c2);
  for(i = 0; i < c1->row; i++){
    for(j = 0; j < c1->col; j++){
      if(FLOAT_EQ(c1->data[i][j], c2->data[i][j], 1e-6) == 1){
        continue;
      }
      else{
        rmsd = 99999.f;
        break;
      }
    }
  }

  DelMatrix(&c1);
  DelMatrix(&c2);
  DelMatrix(&A_T);
  DelMatrix(&A_aligned_T_tmp);
  DelDVector(&t_);
  DelMatrix(&R_);
  return rmsd;
}

double GeomAnalysis::SimplexAlign3DConformations(MOLECULE *m1, MOLECULE *m2)
{
  dvector *steps;
  dvector *x0;
  dvector *best;

  /* 1)Calculate the centroids */
  /* Calculate centroids from molecules
  *  N.B.: molecules must have same number of atoms!
  */
  if(m1->n_atoms != m2->n_atoms){
    printf("[SimplexAlign3DConformations Error]: The number of atoms between conformers differ!\n");
    return -1.;
  }

  NewMatrix(&_Aconf_, m1->n_atoms, 3);
  NewMatrix(&_Bconf_, m2->n_atoms, 3);

  /* 2) Finding the optimal rotation */
  for(int i = 0; i < m1->n_atoms; i++){
    _Aconf_->data[i][0] = m1->atoms[i].coord.x;
    _Aconf_->data[i][1] = m1->atoms[i].coord.y;
    _Aconf_->data[i][2] = m1->atoms[i].coord.z;

    _Bconf_->data[i][0] = m2->atoms[i].coord.x;
    _Bconf_->data[i][1] = m2->atoms[i].coord.y;
    _Bconf_->data[i][2] = m2->atoms[i].coord.z;
  }

  NewDVector(&steps, 6);
  NewDVector(&x0, 6); //ConvertToRandt_v2_
  //NewDVector(&x0, 12);//ConvertToRandt_v1
  DVectorSet(x0, 50.f);
  steps->data[0] = 30.f;
  steps->data[1] = 30.f;
  steps->data[2] = 30.f;
  steps->data[3] = 100.f;
  steps->data[4] = 100.f;
  steps->data[5] = 100.f;
  initDVector(&best);

  double rmsd_min;
  //double rmsd_min = NelderMeadSimplex(&GeomAnalysis::func, x0, NULL, 1e-10, 1000, &best);

  matrix *R;
  dvector *t;
  NewMatrix(&R, 3, 3);
  NewDVector(&t, 3);
  //ConvertToRandt_v1(best, R, t);
  ConvertToRandt_v2_(best, R, t);

  matrix *A_T;
  NewMatrix(&A_T, _Aconf_->col, _Aconf_->row);
  MatrixTranspose(_Aconf_, A_T);

  matrix *A_aligned_T_tmp;
  initMatrix(&A_aligned_T_tmp);

  ComputeCoordinates_(A_T, R, t, &A_aligned_T_tmp);

  DelMatrix(&A_T);

  for(int i = 0; i < m1->n_atoms; i++){
    m1->atoms[i].coord.x = A_aligned_T_tmp->data[0][i];
    m1->atoms[i].coord.y = A_aligned_T_tmp->data[1][i];
    m1->atoms[i].coord.z = A_aligned_T_tmp->data[2][i];
  }

  DelMatrix(&A_aligned_T_tmp);
  /*printf("RMSD: % 3.4f\n", best_rmsd);*/
  DelMatrix(&_Aconf_);
  DelMatrix(&_Bconf_);
  DelMatrix(&R);
  DelDVector(&t);
  DelDVector(&steps);
  return rmsd_min;
}
