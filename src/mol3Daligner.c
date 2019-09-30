/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "mol3Daligner.h"
#include "molecule.h"
#include <scientific.h>
#include <math.h>
#include <sys/timeb.h>


void TranslateConformation(MOLECULE m, double cx, double cy, double cz)
{
  size_t i;
  for(i = 0; i < m.n_atoms; i++){
    m.atoms[i].coord.x = m.atoms[i].coord.x+cx;
    m.atoms[i].coord.y = m.atoms[i].coord.y+cy;
    m.atoms[i].coord.z = m.atoms[i].coord.z+cz;
  }
}


/*
 * pure translation
 */
void mk3dtr(double x, double y, double z, matrix **t)
{
  ResizeMatrix(t, 4, 4);
  (*t)->data[0][0] = (*t)->data[1][1] = (*t)->data[2][2] = (*t)->data[3][3] = 1.f;
  (*t)->data[0][3] = x;
  (*t)->data[1][3] = y;
  (*t)->data[2][3] = z;
}

/*
 * pure rotation around x axis
 */
void mk3drotx(double theta, matrix **xr)
{
  ResizeMatrix(xr, 4, 4);
  (*xr)->data[0][0] = (*xr)->data[1][1] = (*xr)->data[2][2] = (*xr)->data[3][3] = 1.f;
  (*xr)->data[1][1] = cos(theta); (*xr)->data[1][2] = -sin(theta);
  (*xr)->data[2][1] = sin(theta); (*xr)->data[2][2] = cos(theta);
}

/*
 * pure rotation around y axis
 */
void mk3droty(double theta, matrix **yr)
{
  ResizeMatrix(yr, 4, 4);
  (*yr)->data[0][0] = (*yr)->data[1][1] = (*yr)->data[2][2] = (*yr)->data[3][3] = 1.f;
  (*yr)->data[0][0] = cos(theta); (*yr)->data[0][2] = -sin(theta);
  (*yr)->data[2][0] = sin(theta); (*yr)->data[2][2] = cos(theta);
}

/*
 * pure rotation around z axis
 */
void mk3drotz(double theta, matrix **zr)
{
  ResizeMatrix(zr, 4, 4);
  (*zr)->data[0][0] = (*zr)->data[1][1] = (*zr)->data[2][2] = (*zr)->data[3][3] = 1.f;
  (*zr)->data[0][0] = cos(theta); (*zr)->data[0][1] = -sin(theta);
  (*zr)->data[1][0] = sin(theta); (*zr)->data[1][1] = cos(theta);
}

/*
 * Random conformational rotation using quaternions
 */

double randfrom(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void RandomConformationRotation(MOLECULE *m)
{
  size_t i;
  dvector *p, *np;
  matrix *xr, *yr, *zr, *tr, *trb, *zyr, *zyxr;

  initMatrix(&xr);
  initMatrix(&yr);
  initMatrix(&zr);
  initMatrix(&tr);
  initMatrix(&trb);

  double cx, cy, cz;
  cx = cy = cz = 0.f;
  for(i = 0; i < m->n_atoms; i++){
    cx += m->atoms[i].coord.x;
    cy += m->atoms[i].coord.y;
    cz += m->atoms[i].coord.z;
  }

  cx /= (double)m->n_atoms;
  cy /= (double)m->n_atoms;
  cz /= (double)m->n_atoms;

  mk3dtr(-cx, -cy, -cz, &tr);
  mk3dtr(cx/2., cy/2., cz/2., &trb);

  struct timeb t_current;
  ftime(&t_current);
  /*printf("%hu\n", t_current.millitm);*/
  srand(t_current.millitm);

  double yaw, pitch, roll;
  /*yaw = pitch = roll = 0.f;*/

  yaw = randfrom(-_pi_, _pi_);
  pitch = randfrom(-_pi_, _pi_);
  roll = randfrom(-_pi_, _pi_);
  /*printf("%f %f %f\n", yaw, pitch, roll);*/

  mk3drotx(roll, &xr);
  mk3droty(pitch, &yr);
  mk3drotz(yaw, &zr);

  NewMatrix(&zyr, 4, 4);
  NewMatrix(&zyxr, 4, 4);
  /*PrintMatrix(zyxr);*/
  /*Random Molecular Conformation Rotation*/
  NewDVector(&p, 4);
  NewDVector(&np, 4);

  MatrixDotProduct(zr, yr, zyr);
  MatrixDotProduct(zyr, xr, zyxr);

  for(i = 0; i < m->n_atoms; i++){
    p->data[0] = m->atoms[i].coord.x;
    p->data[1] = m->atoms[i].coord.y;
    p->data[2] = m->atoms[i].coord.z;
    p->data[3] = np->data[3] = 1.f;
    /*translate p to np using tr*/
    MatrixDVectorDotProduct(tr, p, np);
    /*reset p*/
    p->data[0] = 0.f;
    p->data[1] = 0.f;
    p->data[2] = 0.f;
    p->data[3] = 0.f;
    /*rotate np to p using zyxr*/
    MatrixDVectorDotProduct(zyxr, np, p);
    /*reset np*/
    np->data[0] = 0.f;
    np->data[1] = 0.f;
    np->data[2] = 0.f;
    np->data[3] = 0.f;
    /*translate back p to np using trb*/
    MatrixDVectorDotProduct(trb, p, np);
    m->atoms[i].coord.x = np->data[0];
    m->atoms[i].coord.y = np->data[1];
    m->atoms[i].coord.z = np->data[2];
    /*reset np*/
    np->data[0] = 0.f;
    np->data[1] = 0.f;
    np->data[2] = 0.f;
    np->data[3] = 0.f;
  }

  DelDVector(&np);
  DelDVector(&p);
  DelMatrix(&xr);
  DelMatrix(&yr);
  DelMatrix(&zr);
  DelMatrix(&zyr);
  DelMatrix(&zyxr);
  DelMatrix(&tr);
  DelMatrix(&trb);
}

/*
* Roto translation work for 2D, 3D, and N-Dimensional problems
* If a dataset A want to be aligned in B there is an equation to solve for R,t:
* B = R*A + t
*
* Where R is a rotation matrix
* and t is a translation vector
*
* Algorithm:
* 1) Find the centroid in both datasets
* 2) Bring both dataset to the origin then find the optimal rotation by
*    using the SVD
* 3) Find the translation t
*/

void RotoTranslation(matrix *A, matrix *B, matrix **R, dvector **t)
{
  dvector *cm1, *cm2;
  matrix *AA_T, *BB, *H; /*used to calculate the covariance matrix H for the step 2 */

  if(A->row != B->row && A->col != B->col){ return; }

  NewMatrix(&AA_T, A->col, A->row);
  MatrixTranspose(A, AA_T);

  NewMatrix(&BB, B->row, B->col);
  MatrixCopy(B, &BB);

  /* 1)Calculate the centroids */

    /* Calculate centroids from molecules
    *  N.B.: molecules must have same number of atoms!
    */

    int n = A->row;
    int m = A->col;
    NewDVector(&cm1, m);
    NewDVector(&cm2, m);

    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        cm1->data[j] += AA_T->data[j][i];
        cm2->data[j] += BB->data[i][j];
      }
    }

    for(int j = 0; j < m; j++){
      cm1->data[j] /= (float)n;
      cm2->data[j] /= (float)n;
    }

    /* 2) Finding the optimal rotation */
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        AA_T->data[j][i] -= cm1->data[j];
        BB->data[i][j] -= cm2->data[j];
      }
    }

    /*Calculate the covariance matrix H */
    NewMatrix(&H, m, m);
    MatrixDotProduct(AA_T, BB, H);

    /*SVD of H */
    matrix *U, *S, *VT;
    initMatrix(&U);
    initMatrix(&S);
    initMatrix(&VT);
    /*SVD(H, &U, &S, &VT); //It's working but not well as lapack version! */
    SVDlapack(H, &U, &S, &VT);

    /*puts("H");
    PrintMatrix(H);
    PrintMatrix(U);
    PrintMatrix(S);
    PrintMatrix(VT);
    puts("-------------");
    */
    matrix *U_T, *V;
    NewMatrix(&U_T, U->col, U->row);
    MatrixTranspose(U, U_T);
    NewMatrix(&V, VT->col, VT->row);
    MatrixTranspose(VT, V);

    /* Calculate the rotation R */
    ResizeMatrix(R, m, m);
    MatrixDotProduct(V, U_T, (*R));

    /* Now R is the rotation matrix! */
    /* Check for the special reflection case */

    if(MatrixDeterminant((*R)) < 0){
       /*Reflection detected! */
      for(int i = 0; i < V->row; i++){
        V->data[i][m-1] *= -1;
      }

      MatrixSet((*R), 0.f);
      MatrixDotProduct(V, U_T, (*R));
    }

    DelMatrix(&U_T);
    DelMatrix(&V);
    DelMatrix(&S);
    DelMatrix(&U);
    DelMatrix(&VT);
    DelMatrix(&H);

    /* 3) Calculate translation
    t = -R x centroid A  + centroid B
    */

    matrix *mR;
    NewMatrix(&mR, m, m);
    for(int i = 0; i < (*R)->row; i++){
      for(int j = 0; j < (*R)->col; j++)
        mR->data[i][j] = (*R)->data[i][j] *-1;
    }

    DVectorResize(t, m);
    MatrixDVectorDotProduct(mR, cm1, (*t));
    for(int i = 0; i < (*t)->size; i++)
      (*t)->data[i] += cm2->data[i];

    DelMatrix(&mR);
    DelMatrix(&AA_T);
    DelMatrix(&BB);
    DelDVector(&cm1);
    DelDVector(&cm2);

    PrintMatrix((*R));
    PrintDVector((*t));
}

double RMSD(matrix *A, matrix *B)
{
  size_t i, j;
  double rmsd = 0.f;
  if(A->col == B->col && A->row == B->row){
    for(j = 0; j < A->col; j++){
      for(i = 0; i < A->row; i++){
        rmsd += square(A->data[i][j] - B->data[i][j]);
      }
    }
    rmsd = sqrt(rmsd/(float)A->row);
  }
  else if(A->col == B->row && A->row == B->col){ // if a or b are tranposed
    for(j = 0; j < A->col; j++){
      for(i = 0; i < A->row; i++){
        rmsd += square(A->data[i][j] - B->data[j][i]);
      }
    }
    rmsd = sqrt(rmsd/(float)A->col);
  }
  else{
    rmsd = 9999.f;
  }
  return rmsd;
}

void ComputeCoordinates(matrix *A_T, matrix *R, dvector *t, matrix **A_aligned_T)
{
  ResizeMatrix(A_aligned_T, A_T->row, A_T->col);

  MatrixDotProduct(R, A_T, (*A_aligned_T));

  for(int i = 0; i < A_T->col; i++){
    for(int j = 0; j < A_T->row; j++)
      (*A_aligned_T)->data[j][i] += t->data[j];
  }
}

double Align3DPharmacophore(MOLECULE m1, MOLECULE m2, uivector *aid1, uivector *aid2)
{
  matrix *A, *B, *R, *A_aligned_T;
  dvector *t;
  double best_rmsd = 9999.;
  matrix *R_best;
  dvector *t_best;
  int best_it = 0;

  /* To align the whole m1 molecule we need the rotation matrix R and translation vector t */
  initMatrix(&R_best);
  initDVector(&t_best);

  /* 1)Calculate the centroids */
  /* Calculate centroids from molecules
  *  N.B.: molecules must have same number of atoms!
  */
  if(aid1->size != aid2->size){
    printf("[Align3DPharmacophore Error]: The number of pharmacophoric atoms between molecules differ!\n");
    return -9999.;
  }

  NewMatrix(&A, aid1->size, 3);
  NewMatrix(&B, aid2->size, 3);

  /* 2) Finding the optimal rotation */
  /* Select the pharmacophoric atoms to align */
  for(int i = 0; i < aid1->size; i++){
    A->data[i][0] = m1.atoms[aid1->data[i]-1].coord.x;
    A->data[i][1] = m1.atoms[aid1->data[i]-1].coord.y;
    A->data[i][2] = m1.atoms[aid1->data[i]-1].coord.z;

    B->data[i][0] = m2.atoms[aid2->data[i]-1].coord.x;
    B->data[i][1] = m2.atoms[aid2->data[i]-1].coord.y;
    B->data[i][2] = m2.atoms[aid2->data[i]-1].coord.z;

    printf("%s with %s\n",  m1.atoms[aid1->data[i]-1].type, m2.atoms[aid2->data[i]-1].type);
  }


  matrix *A_p;
  NewMatrix(&A_p, A->row, A->col);
  /* run the calcuation doing different row permutation in A_T */
  int p[6][3] = {{1,2,3}, {1,3,2}, {2,1,3}, {2,3,1}, {3,1,2}, {3,2,1}};

  int iter = 0;
  while(iter < 100){
    for(int it = 0; it < 6; it++){
      initMatrix(&R);
      initDVector(&t);
      for(int i = 0; i < A->row; i++){
        for(int j = 0; j < A->col; j++){
          A_p->data[i][j] = A->data[i][p[it][j]-1];
        }
      }

      /* fit A to B so calculate the rotation matrix R and the translation vector t*/
      RotoTranslation(A_p, B, &R, &t);

      /* Now change the coordinates of the pharmacophore 1 using R and t
       * in order to make effective the alignment to the pharmacophore 2*/
      matrix *A_T;
      NewMatrix(&A_T, A->col, A->row);
      MatrixTranspose(A_p, A_T);

      matrix *A_aligned_T_tmp;
      initMatrix(&A_aligned_T_tmp);

      ComputeCoordinates(A_T, R, t, &A_aligned_T_tmp);

      DelMatrix(&A_T);

      double rmsd = RMSD(A_aligned_T_tmp, B);

      /*printf("TMP RMSD iter %d %f\n", iter, rmsd);*/
      if(rmsd < best_rmsd){
         best_rmsd = rmsd;
         MatrixCopy(R, &R_best);
         DVectorCopy(t, &t_best);
         best_it = it;
      }
      DelMatrix(&A_aligned_T_tmp);
      DelDVector(&t);
      DelMatrix(&R);
    }
    iter++;
  }

  DelMatrix(&A_p);
  DelMatrix(&A);
  DelMatrix(&B);

  /* Now change the coordinates to m1 with R and t
   * in order to make effective the alignment to m2 */
  matrix *Am1;
  NewMatrix(&Am1, m1.n_atoms, 3);

  for(int i = 0; i < m1.n_atoms; i++){
    Am1->data[i][0] = m1.atoms[i].coord.x;
    Am1->data[i][1] = m1.atoms[i].coord.y;
    Am1->data[i][2] = m1.atoms[i].coord.z;
  }

  matrix *A_T;
  NewMatrix(&A_T, 3, m1.n_atoms);

  for(int i = 0; i < m1.n_atoms; i++){
    for(int j = 0; j < 3; j++){
      A_T->data[j][i] = Am1->data[i][p[best_it][j]-1];
    }
  }
  DelMatrix(&Am1);

  NewMatrix(&A_aligned_T, 3, m1.n_atoms);

  ComputeCoordinates(A_T, R_best, t_best, &A_aligned_T);

  DelMatrix(&A_T);

  for(int i = 0; i < m1.n_atoms; i++){
    m1.atoms[i].coord.x = A_aligned_T->data[0][i];
    m1.atoms[i].coord.y = A_aligned_T->data[1][i];
    m1.atoms[i].coord.z = A_aligned_T->data[2][i];
  }

  /*printf("RMSD: % 3.4f\n", best_rmsd);*/
  DelMatrix(&A_aligned_T);
  DelMatrix(&R_best);
  DelDVector(&t_best);

  return best_rmsd;
}

double Align3DConformations(MOLECULE m1, MOLECULE m2)
{
  matrix *A, *B, *R, *A_aligned_T;
  dvector *t;
  double best_rmsd = 9999.;


  NewMatrix(&A_aligned_T, 3, m1.n_atoms);

  /* 1)Calculate the centroids */
  /* Calculate centroids from molecules
  *  N.B.: molecules must have same number of atoms!
  */
  if(m1.n_atoms != m2.n_atoms){
    printf("[Align3DConformations Error]: The number of atoms between conformers differ!\n");
    return -1.;
  }

  NewMatrix(&A, m1.n_atoms, 3);
  NewMatrix(&B, m2.n_atoms, 3);

  /* 2) Finding the optimal rotation */
  for(int i = 0; i < m1.n_atoms; i++){
    A->data[i][0] = m1.atoms[i].coord.x;
    A->data[i][1] = m1.atoms[i].coord.y;
    A->data[i][2] = m1.atoms[i].coord.z;

    B->data[i][0] = m2.atoms[i].coord.x;
    B->data[i][1] = m2.atoms[i].coord.y;
    B->data[i][2] = m2.atoms[i].coord.z;
  }

  matrix *A_p;
  NewMatrix(&A_p, A->row, A->col);
  /* run the calcuation doing different row permutation in A_T */
  int p[6][3] = {{1,2,3}, {1,3,2}, {2,1,3}, {2,3,1}, {3,1,2}, {3,2,1}};

  for(int it = 0; it < 6; it++){
    initMatrix(&R);
    initDVector(&t);
    for(int i = 0; i < A->row; i++){
      for(int j = 0; j < A->col; j++){
        A_p->data[i][j] = A->data[i][p[it][j]-1];
      }
    }

    /* fit A to B so calculate the rotation matrix R and the translation vector t*/
    RotoTranslation(A_p, B, &R, &t);

    /* Now change the coordinates to m1 with R and t
     * in order to make effective the alignment to m2 */
    matrix *A_T;
    NewMatrix(&A_T, A->col, A->row);
    MatrixTranspose(A_p, A_T);

    matrix *A_aligned_T_tmp;
    initMatrix(&A_aligned_T_tmp);

    ComputeCoordinates(A_T, R, t, &A_aligned_T_tmp);

    DelMatrix(&A_T);

    double rmsd = RMSD(A_aligned_T_tmp, B);
    //printf("TMP RMSD %f\n", rmsd);
    if(rmsd < best_rmsd){
       best_rmsd = rmsd;
       MatrixCopy(A_aligned_T_tmp, &A_aligned_T);
    }
    DelMatrix(&A_aligned_T_tmp);
    DelDVector(&t);
    DelMatrix(&R);
  }

  DelMatrix(&A_p);
  DelMatrix(&A);
  DelMatrix(&B);


  for(int i = 0; i < m1.n_atoms; i++){
    m1.atoms[i].coord.x = A_aligned_T->data[0][i];
    m1.atoms[i].coord.y = A_aligned_T->data[1][i];
    m1.atoms[i].coord.z = A_aligned_T->data[2][i];
  }

  /*printf("RMSD: % 3.4f\n", best_rmsd);*/
  DelMatrix(&A_aligned_T);
  return best_rmsd;
}

void ConvertToRandt_v1(dvector *x, matrix *R, dvector *t)
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

void ConvertToRandt_v2(dvector *x, matrix *R, dvector *t)
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

double func(dvector *x)
{
  size_t i, j;
  matrix *R_;
  dvector *t_;
  NewMatrix(&R_, 3, 3);
  NewDVector(&t_, 3);

  PrintDVector(x);
  //ConvertToRandt_v1(x, R_, t_);
  ConvertToRandt_v2(x, R_, t_);

  /* Now change the coordinates to m1 with R and t
   * in order to make effective the alignment to m2 */
  matrix *A_T;
  NewMatrix(&A_T, _Aconf_->col, _Aconf_->row);
  MatrixTranspose(_Aconf_, A_T);

  matrix *A_aligned_T_tmp;
  initMatrix(&A_aligned_T_tmp);

  ComputeCoordinates(A_T, R_, t_, &A_aligned_T_tmp);

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

/*~Not working well!!! molecular deformation */
double SimplexAlign3DConformations(MOLECULE m1, MOLECULE m2)
{
  dvector *steps;
  dvector *x0;
  dvector *best;

  /* 1)Calculate the centroids */
  /* Calculate centroids from molecules
  *  N.B.: molecules must have same number of atoms!
  */
  if(m1.n_atoms != m2.n_atoms){
    printf("[SimplexAlign3DConformations Error]: The number of atoms between conformers differ!\n");
    return -1.;
  }

  NewMatrix(&_Aconf_, m1.n_atoms, 3);
  NewMatrix(&_Bconf_, m2.n_atoms, 3);

  /* 2) Finding the optimal rotation */
  for(int i = 0; i < m1.n_atoms; i++){
    _Aconf_->data[i][0] = m1.atoms[i].coord.x;
    _Aconf_->data[i][1] = m1.atoms[i].coord.y;
    _Aconf_->data[i][2] = m1.atoms[i].coord.z;

    _Bconf_->data[i][0] = m2.atoms[i].coord.x;
    _Bconf_->data[i][1] = m2.atoms[i].coord.y;
    _Bconf_->data[i][2] = m2.atoms[i].coord.z;
  }

  NewDVector(&steps, 6);
  NewDVector(&x0, 6); //ConvertToRandt_v2
  //NewDVector(&x0, 12);//ConvertToRandt_v1
  DVectorSet(x0, 50.f);
  steps->data[0] = 30.f;
  steps->data[1] = 30.f;
  steps->data[2] = 30.f;
  steps->data[3] = 100.f;
  steps->data[4] = 100.f;
  steps->data[5] = 100.f;
  initDVector(&best);

  double rmsd_min = NelderMeadSimplex(&func, x0, NULL, 1e-10, 1000, &best);

  matrix *R;
  dvector *t;
  NewMatrix(&R, 3, 3);
  NewDVector(&t, 3);
  //ConvertToRandt_v1(best, R, t);
  ConvertToRandt_v2(best, R, t);

  matrix *A_T;
  NewMatrix(&A_T, _Aconf_->col, _Aconf_->row);
  MatrixTranspose(_Aconf_, A_T);

  matrix *A_aligned_T_tmp;
  initMatrix(&A_aligned_T_tmp);

  ComputeCoordinates(A_T, R, t, &A_aligned_T_tmp);

  DelMatrix(&A_T);

  for(int i = 0; i < m1.n_atoms; i++){
    m1.atoms[i].coord.x = A_aligned_T_tmp->data[0][i];
    m1.atoms[i].coord.y = A_aligned_T_tmp->data[1][i];
    m1.atoms[i].coord.z = A_aligned_T_tmp->data[2][i];
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


double ComputeRMSD(MOLECULE m1, MOLECULE m2)
{
  int i;
  matrix *A, *B;
  double rmsd = 9999.f;
  NewMatrix(&A, m1.n_atoms, 3);
  NewMatrix(&B, m2.n_atoms, 3);

  /* 2) Finding the optimal rotation */
  for(i = 0; i < m1.n_atoms; i++){
    A->data[i][0] = m1.atoms[i].coord.x;
    A->data[i][1] = m1.atoms[i].coord.y;
    A->data[i][2] = m1.atoms[i].coord.z;

    B->data[i][0] = m2.atoms[i].coord.x;
    B->data[i][1] = m2.atoms[i].coord.y;
    B->data[i][2] = m2.atoms[i].coord.z;
  }

  rmsd = RMSD(A, B);
  DelMatrix(&A);
  DelMatrix(&B);
  return rmsd;
}
