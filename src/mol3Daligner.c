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
#include "shapedesc.h"
#include <scientific.h>
#include <math.h>
#include <sys/time.h>

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

  struct timeval tv;
  gettimeofday(&tv,NULL);
  /*printf("%hu\n", t_current.millitm);*/
  srand(tv.tv_usec);

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

void Random3DMatrixRotation(matrix *A)
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
  for(i = 0; i < A->row; i++){
    cx += A->data[i][0];
    cy += A->data[i][1];
    cz += A->data[i][2];
  }

  cx /= (double)A->row;
  cy /= (double)A->row;
  cz /= (double)A->row;

  mk3dtr(-cx, -cy, -cz, &tr);
  mk3dtr(cx/2., cy/2., cz/2., &trb);

  struct timeval tv;
  gettimeofday(&tv,NULL);
  /*printf("%hu\n", t_current.millitm);*/
  srand(tv.tv_usec);

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

  for(i = 0; i < A->row; i++){
    p->data[0] = A->data[i][0];
    p->data[1] = A->data[i][1];
    p->data[2] = A->data[i][2];
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
    A->data[i][0] = np->data[0];
    A->data[i][1] = np->data[1];
    A->data[i][2] = np->data[2];
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

void find_nn_points(matrix *src_pts,
                    matrix *dst_pts,
                    matrix **nn_pts)
{
  /* K-NN algorithm style where K is 1
   * This algorithm give back the index of dst that are nearest to src
   */
  size_t i, j;
  matrix *dstmx;
  initMatrix(&dstmx);
  EuclideanDistance(dst_pts, src_pts, &dstmx, 1);

  /* dstmx
   * Every row in dstmx correspond to the src->row
   * every column in dstmx correspond to dst->row
   */
  uivector *nn_indx;
  NewUIVector(&nn_indx, dstmx->row);

  for(i = 0; i < dstmx->row; i++){
    int min_j = 0;
    for(int j = 1; j < dstmx->col; j++){
      if(dstmx->data[i][j] < dstmx->data[i][min_j]){
        min_j = j;
      }
    }
    nn_indx->data[i] = min_j;
  }


  for(i = 0; i < nn_indx->size; i++){
    for(j = 0; j < dst_pts->col; j++){
      (*nn_pts)->data[i][j] = dst_pts->data[nn_indx->data[i]][j];
    }
  }

  DelUIVector(&nn_indx);
  DelMatrix(&dstmx);
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
    /*SVD(H, &U, &S, &VT);*/
    SVDlapack(H, &U, &S, &VT);

    /*puts("H");
    PrintMatrix(H);
    PrintMatrix(U);
    PrintMatrix(S);
    PrintMatrix(VT);

    puts("-------------");
    abort();
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
      for(int j = 0; j < VT->col; j++){
        VT->data[VT->row-1][j] *=-1;
      }
      MatrixSet(V, 0.f);
      MatrixTranspose(VT, V);
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
    t = centroid B.T  - (R* centroid A.T)
    */
    dvector *Rcm;
    NewDVector(&Rcm, m);
    MatrixDVectorDotProduct((*R), cm1, Rcm);
    DVectorResize(t, m);
    for(size_t i = 0; i < (*t)->size; i++){
      (*t)->data[i] = cm2->data[i] - Rcm->data[i];
    }
    DelDVector(&Rcm);

    DelMatrix(&AA_T);
    DelMatrix(&BB);
    DelDVector(&cm1);
    DelDVector(&cm2);

    /*
    PrintMatrix((*R));
    PrintDVector((*t));
    */
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

void ApplyTransformation(matrix *pts,
                         matrix *R,
                         dvector *t,
                         matrix **pts_transformed)
{
  size_t i, j;
  matrix *pts_T, *pts_transformed_T;
  NewMatrix(&pts_T, pts->col, pts->row);
  NewMatrix(&pts_transformed_T, pts->col, pts->row);
  MatrixTranspose(pts, pts_T);
  ResizeMatrix(pts_transformed, pts->row, pts->col);
  MatrixDotProduct(R, pts_T, pts_transformed_T);
  for(i = 0; i < pts->row; i++){
    for(j = 0; j < pts->col; j++){
      pts_transformed_T->data[j][i] += t->data[j];
    }
  }

  MatrixTranspose(pts_transformed_T, (*pts_transformed));
  DelMatrix(&pts_transformed_T);
  DelMatrix(&pts_T);
}

void RotoTranslateMolecule(MOLECULE m, matrix *R, dvector *t)
{
  matrix *A;
  matrix *A_aligned_T;
  NewMatrix(&A, m.n_atoms, 3);

  for(int i = 0; i < m.n_atoms; i++){
    A->data[i][0] = m.atoms[i].coord.x;
    A->data[i][1] = m.atoms[i].coord.y;
    A->data[i][2] = m.atoms[i].coord.z;
  }

  matrix *A_T;
  NewMatrix(&A_T, 3, m.n_atoms);
  MatrixTranspose(A, A_T);
  DelMatrix(&A);

  NewMatrix(&A_aligned_T, 3, m.n_atoms);
  ComputeCoordinates(A_T, R, t, &A_aligned_T);
  DelMatrix(&A_T);

  for(int i = 0; i < m.n_atoms; i++){
    m.atoms[i].coord.x = A_aligned_T->data[0][i];
    m.atoms[i].coord.y = A_aligned_T->data[1][i];
    m.atoms[i].coord.z = A_aligned_T->data[2][i];
  }
  DelMatrix(&A_aligned_T);
}

double ICP(matrix *src_pts_,
           matrix *dst_pts,
           double tolerance,
           matrix **R,
           dvector **t,
           MOLECULE *src_mol)
{
    matrix *R_;
    dvector *t_;
    matrix *src_pts;
    matrix *current_pts;
    matrix *neigh_pts;
    NewMatrix(&src_pts, src_pts_->row, src_pts_->col);
    MatrixCopy(src_pts_, &src_pts);

    NewMatrix(&current_pts, src_pts->row, src_pts->col);
    NewMatrix(&neigh_pts, src_pts->row, src_pts->col);
    MatrixCopy(src_pts, &current_pts);


    double last_rmsd = 9999.f;
    size_t k = 0;
    while(1){
      find_nn_points(current_pts, dst_pts, &neigh_pts);
      initMatrix(&R_);
      initDVector(&t_);
      RotoTranslation(src_pts, neigh_pts, &R_, &t_);
      ApplyTransformation(src_pts, R_, t_, &current_pts);
      if(src_mol != NULL){
        /*debugging purpose*/
        char buf[24];
        snprintf(buf, 24, "tmp_%zu_mol_aligned.mol2", k);
        for(int i = 0; i < (*src_mol).n_atoms; i++){
          (*src_mol).atoms[i].coord.x = current_pts->data[i][0];
          (*src_mol).atoms[i].coord.y = current_pts->data[i][1];
          (*src_mol).atoms[i].coord.z = current_pts->data[i][2];
        }        
        SaveMol2Molecule((*src_mol), buf);
      }
      double rmsd = RMSD(current_pts, neigh_pts);
      if(fabs(rmsd-last_rmsd) < tolerance){
        MatrixCopy(R_, R);
        DVectorCopy(t_, t);
        DelMatrix(&R_);
        DelDVector(&t_);
        break;
      }
      else{
        //printf("%f -> %f\n", last_rmsd, rmsd);
        last_rmsd = rmsd;
        MatrixCopy(R_, R);
        DVectorCopy(t_, t);
        DelMatrix(&R_);
        DelDVector(&t_);
      }
      k +=1;
    }

    DelMatrix(&src_pts);
    DelMatrix(&neigh_pts);
    DelMatrix(&current_pts);
    return last_rmsd;
}


/* get3DAnchorPoints act on 3D coordinates */
void get3DAnchorPoints(matrix *c, uivector *aidx, size_t npnts, uivector **aid)
{
 /*
  * c is the matrix with xyz coordinates
  * aidx is a map of atom id in c and the real id on the molecule
  * aid is the vector of point to be filled up
  */

  size_t i, j, k;
  matrix *points;
  /*make copy of c and aidx*/
  matrix *cc;
  NewMatrix(&cc, c->row, c->col);
  MatrixCopy(c, &cc);
  uivector *caidx;
  NewUIVector(&caidx, aidx->size);
  for(i = 0; i < aidx->size; i++){
    caidx->data[i] = aidx->data[i];
  }
  /*Calcualte baricenter*/
  NewMatrix(&points, 1, 3);
  for(i = 0; i < c->row; i++){
    points->data[0][0] += c->data[i][0];
    points->data[0][1] += c->data[i][1];
    points->data[0][2] += c->data[i][2];
  }
  points->data[0][0] /= (double)c->row;
  points->data[0][1] /= (double)c->row;
  points->data[0][2] /= (double)c->row;

  /* Select the points */
  for(i = 0; i < npnts; i++){
    matrix *dst;
    initMatrix(&dst);
    EuclideanDistance(cc, points, &dst, 1);
    dvector *dst_sum;
    NewDVector(&dst_sum, c->row);
    for(j = 0; j < dst->row; j++){
      for(k = 0; k < dst->col; k++){
        dst_sum->data[k] += dst->data[j][k];
      }
    }

    /*Get the point */
    double max_dst = dst_sum->data[j];
    size_t max_dst_indx = 0;
    for(j = 1; j < dst_sum->size; j++){
      if(dst_sum->data[j] > max_dst){
        max_dst = dst_sum->data[j];
        max_dst_indx = j;
      }
      else{
        continue;
      }
    }

    UIVectorAppend(aid, caidx->data[max_dst_indx]);
    dvector *newrow;
    NewDVector(&newrow, 3);
    newrow->data[0] = c->data[max_dst_indx][0];
    newrow->data[1] = c->data[max_dst_indx][1];
    newrow->data[2] = c->data[max_dst_indx][2];
    MatrixAppendRow(&points, newrow);
    DelDVector(&newrow);
    DelMatrix(&dst);
    DelDVector(&dst_sum);
    UIVectorRemoveAt(&caidx, max_dst_indx);
    MatrixDeleteRowAt(&cc, max_dst_indx);
  }
  DelMatrix(&points);
  DelMatrix(&cc);
  DelUIVector(&caidx);
}

void get3DAnchorPointsMaxDis(matrix *c, uivector *aidx, size_t npoints, uivector **aid)
{
 /*
  * c is the matrix with xyz coordinates
  * aidx is a map of atom id in c and the real id on the molecule
  * aid is the vector of point to be filled up
  */

  int run = SIGSCIENTIFICRUN;
  uivector *pnt;
  initUIVector(&pnt);
  MaxDis(c,  npoints, 0, &pnt, 4, &run);

  UIVectorResize(aid, npoints);
  for(int i = 0; i < pnt->size; i++){
    (*aid)->data[i] = aidx->data[pnt->data[i]];
  }
  DelUIVector(&pnt);
}

/* getAnchorPoints act on 2D pca plane*/
void getAnchorPoints(matrix *scores, uivector *aidx, uivector **aid)
{
  double minx, miny, maxx, maxy;
  minx = scores->data[0][0];
  maxx = scores->data[0][0];
  miny = scores->data[0][1];
  maxy = scores->data[0][1];

  UIVectorResize(aid, 4);
  (*aid)->data[0] = aidx->data[0];
  (*aid)->data[1] = aidx->data[0];
  (*aid)->data[2] = aidx->data[0];
  (*aid)->data[3] = aidx->data[0];
  for(int i = 1; i < scores->row; i++){
    if(scores->data[i][0] < minx){
      minx = scores->data[i][0];
      (*aid)->data[0] = aidx->data[i];
    }

    if(scores->data[i][0] > maxx){
      maxx = scores->data[i][0];
      (*aid)->data[1] = aidx->data[i];
    }

    if(scores->data[i][1] < miny){
      miny = scores->data[i][0];
      (*aid)->data[2] = aidx->data[i];
    }

    if(scores->data[i][1] > maxy){
      maxy = scores->data[i][0];
      (*aid)->data[3] = aidx->data[i];
    }
  }
}

/*DEPRECATED*/
void getMatrixof_H_Coordinates(MOLECULE m, matrix **c, uivector **aidx)
{
  int row = 0;
  for(int i = 0; i <  m.n_atoms; i++){
    if(strcmp( m.atoms[i].asymbl, "H") == 0){
      row++;
    }
    else{
      continue;
    }
  }
  ResizeMatrix(c, row, 3);
  UIVectorResize(aidx, row);
  row = 0;
  for(int i = 0; i < m.n_atoms; i++){
    if(strcmp(m.atoms[i].type, "H") == 0){
      (*c)->data[row][0] = m.atoms[i].coord.x;
      (*c)->data[row][1] = m.atoms[i].coord.y;
      (*c)->data[row][2] = m.atoms[i].coord.z;
      (*aidx)->data[row] = i+1;
      row++;
    }
    else{
      continue;
    }
  }
}
/*DEPRECATED*/
void getMatrixof_Car_Coordinates(MOLECULE m, matrix **c, uivector **aidx)
{
  int row = 0;
  for(int i = 0; i <  m.n_atoms; i++){
    if(strcmp( m.atoms[i].type, "C.ar") == 0){
      row++;
    }
    else{
      continue;
    }
  }
  ResizeMatrix(c, row, 3);
  UIVectorResize(aidx, row);
  row = 0;
  for(int i = 0; i < m.n_atoms; i++){
    if(strcmp(m.atoms[i].type, "C.ar") == 0){
      (*c)->data[row][0] = m.atoms[i].coord.x;
      (*c)->data[row][1] = m.atoms[i].coord.y;
      (*c)->data[row][2] = m.atoms[i].coord.z;
      (*aidx)->data[row] = i+1;
      row++;
    }
    else{
      continue;
    }
  }
}


void getMatrixofAtomTypeChargeScaled(MOLECULE m,
                                     char *atype,
                                     matrix **c,
                                     uivector **aidx)
{
  int row = 0;
  for(int i = 0; i <  m.n_atoms; i++){
    if(strcmp( m.atoms[i].type, atype) == 0){
      row++;
    }
    else{
      continue;
    }
  }
  ResizeMatrix(c, row, 3);
  UIVectorResize(aidx, row);
  row = 0;
  double a = 1.6021765314e-29;
  double b = 3.33564095198152e-30;
  /* We convert these value in  Coulomb * metro (C * m)
   * e Angstrom * 1.6021765314E-19 (C) * 1.0E-10(A) = C*m
   * C*m is converted in debye dividing by 3.336·10^(-30) C·m
   * because 1D = 3.336·10^(-30) C·m
   */
  for(int i = 0; i < m.n_atoms; i++){
    if(strcmp(m.atoms[i].type, atype) == 0){
      (*c)->data[row][0] = (m.atoms[i].coord.x*m.atoms[i].charge*a)/b;
      (*c)->data[row][1] = (m.atoms[i].coord.y*m.atoms[i].charge*a)/b;
      (*c)->data[row][2] = (m.atoms[i].coord.z*m.atoms[i].charge*a)/b;
      (*aidx)->data[row] = i+1;
      row++;
    }
    else{
      continue;
    }
  }
}

double Align3DShapes(MOLECULE m1, MOLECULE m2)
{

  double final_rmse;
  matrix *c1;
  matrix *c2;
  uivector *aid1, *aidx1;
  uivector *aid2, *aidx2;

  initUIVector(&aid1);
  initUIVector(&aid2);

  /* Select 6 shape points */
  initMatrix(&c1);
  initMatrix(&c2);
  initUIVector(&aidx1);
  initUIVector(&aidx2);

  getMatrixof_H_Coordinates(m1, &c1, &aidx1);
  getMatrixof_H_Coordinates(m2, &c2, &aidx2);

  get3DAnchorPoints(c1, aidx1, 6, &aid1);
  get3DAnchorPoints(c2, aidx2, 6, &aid2);

  DelMatrix(&c1);
  DelMatrix(&c2);
  DelUIVector(&aidx1);
  DelUIVector(&aidx2);

  /* Select 4 aromatic points */
  initMatrix(&c1);
  initMatrix(&c2);
  initUIVector(&aidx1);
  initUIVector(&aidx2);

  getMatrixof_Car_Coordinates(m1, &c1, &aidx1);
  getMatrixof_Car_Coordinates(m2, &c2, &aidx2);

  get3DAnchorPoints(c1, aidx1, 4, &aid1);
  get3DAnchorPoints(c2, aidx2, 4, &aid2);

  DelMatrix(&c1);
  DelMatrix(&c2);
  DelUIVector(&aidx1);
  DelUIVector(&aidx2);


  /*Final point selection*/
  //PrintUIVector(aid1);
  //PrintUIVector(aid2);
  /* Once we get the 4 points we select the best
    * alignment by flipping them on the plane
    */
  dvector *rmse;
  NewDVector(&rmse, 4);
  size_t tmpid;

  /* a0 b0
    * a1 b1
    * a2 b2
    * a3 b3
    * No Flip
    */
  rmse->data[0] = Align3DPharmacophore(m1, m2, aid1, aid2);

  /* a0 b1
    * a1 b0
    * a2 b2
    * a3 b3
    * Flip X coordinates
    */
  tmpid = aid2->data[0];
  aid2->data[0] = aid2->data[1];
  aid2->data[1] = tmpid;
  rmse->data[1] = Align3DPharmacophore(m1, m2, aid1, aid2);

  /* a0 b1
    * a1 b0
    * a2 b3
    * a3 b2
    * Flip Y coordinates
    */
  tmpid = aid2->data[2];
  aid2->data[2] = aid2->data[3];
  aid2->data[3] = tmpid;
  rmse->data[2] = Align3DPharmacophore(m1, m2, aid1, aid2);

  /* a0 b0
    * a1 b1
    * a2 b3
    * a3 b2
    * Revert X coordinates
    */
  tmpid = aid2->data[0];
  aid2->data[0] = aid2->data[1];
  aid2->data[1] = tmpid;
  rmse->data[3] = Align3DPharmacophore(m1, m2, aid1, aid2);

  //PrintDVector(rmse);

  int min_rmse = 0;
  for(int i = 1; i < rmse->size; i++){
    if(rmse->data[i] < min_rmse){
      min_rmse = i;
    }
    else{
      continue;
    }
  }

  /*From the last status now... */
  if(min_rmse == 0){
    /*Revert Y*/
    tmpid = aid2->data[2];
    aid2->data[2] = aid2->data[3];
    aid2->data[3] = tmpid;
  }
  else if(min_rmse == 1){
    /*Revert Y and Flip X*/
    tmpid = aid2->data[2];
    aid2->data[2] = aid2->data[3];
    aid2->data[3] = tmpid;

    tmpid = aid2->data[0];
    aid2->data[0] = aid2->data[1];
    aid2->data[1] = tmpid;
  }
  else if(min_rmse == 2){
    /* Flip X */
    tmpid = aid2->data[0];
    aid2->data[0] = aid2->data[1];
    aid2->data[1] = tmpid;
  }
  /*
  else{
      As is...
  }
  */

  final_rmse = Align3DPharmacophore(m1, m2, aid1, aid2);

  DelUIVector(&aid1);
  DelUIVector(&aid2);
  //DelPCAModel(&mod1);
  //DelPCAModel(&mod2);
  DelDVector(&rmse);
  return final_rmse;
}

double Align3DOnVDWShapes(MOLECULE m1, MOLECULE m2)
{
  matrix *A, *B, *R;
  dvector *t;
  double rmsd, last_rmsd =9999.f;
  SHAPEPNT *shap1;
  SHAPEPNT *shap2;

  getShapePoints(m1, &shap1, vanderwaals);
  getShapePoints(m2, &shap2, vanderwaals);

  //WriteMol2SurfPoints(shap1, "mol1.mol2");
  //WriteMol2SurfPoints(shap2, "mol2.mol2");

  uivector *aid1, *aidx1;
  uivector *aid2, *aidx2;

  NewUIVector(&aidx1, shap1->surf->row);
  NewUIVector(&aidx2, shap2->surf->row);
  for(int i = 0; i < shap1->surf->row; i++)
    aidx1->data[i] = i;
  for(int i = 0; i < shap2->surf->row; i++)
    aidx2->data[i] = i;

  int ntpnts = 100;
  //PrintMatrix(shap1->surf);
  //PrintMatrix(shap2->surf);

  //Center to zero the shape points
  MeanCenteredMatrix(shap1->surf, shap1->surf);
  MeanCenteredMatrix(shap2->surf, shap2->surf);

  initUIVector(&aid1);
  initUIVector(&aid2);

  /*sampling the shape*/
  get3DAnchorPointsMaxDis(shap1->surf, aidx1, ntpnts, &aid1);
  get3DAnchorPointsMaxDis(shap2->surf, aidx2, ntpnts, &aid2);

  DelUIVector(&aidx1);
  DelUIVector(&aidx2);


  /* 1)Calculate the centroids */
  /* Calculate centroids from molecules
  *  N.B.: molecules must have same number of atoms!
  */

  NewMatrix(&A, aid1->size, 3);
  NewMatrix(&B, aid2->size, 3);

  /* 2) Finding the optimal rotation */
  /* Select the pharmacophoric atoms to align */
  for(int i = 0; i < aid1->size; i++){
    A->data[i][0] = shap1->surf->data[aid1->data[i]][0];
    A->data[i][1] = shap1->surf->data[aid1->data[i]][1];
    A->data[i][2] = shap1->surf->data[aid1->data[i]][2];

    B->data[i][0] = shap2->surf->data[aid2->data[i]][0];
    B->data[i][1] = shap2->surf->data[aid2->data[i]][1];
    B->data[i][2] = shap2->surf->data[aid2->data[i]][2];
  }

  DelUIVector(&aid1);
  DelUIVector(&aid2);

  while(1){
    initMatrix(&R);
    initDVector(&t);
    rmsd = ICP(A, B, 1e-6, &R, &t, NULL);
    if(fabs(rmsd-last_rmsd) < 1e-6){
      RotoTranslateMolecule(m1, R, t);
      break;
    }
    if(rmsd < last_rmsd)
      last_rmsd = rmsd;
    Random3DMatrixRotation(A);
    DelMatrix(&R);
    DelDVector(&t);
  }

  DelMatrix(&A);
  DelMatrix(&B);

  deleteShapePoints(&shap1);
  deleteShapePoints(&shap2);
  return rmsd;
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

    //printf("%s with %s\n",  m1.atoms[aid1->data[i]-1].type, m2.atoms[aid2->data[i]-1].type);
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
  matrix *A, *B, *R, *A_best;
  dvector *t;
  double rmsd, best_rmsd = 9999.f;

  if(m1.n_atoms != m2.n_atoms){
    printf("[Align3DConformations Error]: The number of atoms between conformers differ!\n");
    return -1.;
  }
  
  NewMatrix(&A, m1.n_atoms, 3);
  NewMatrix(&B, m2.n_atoms, 3);
  NewMatrix(&A_best, m1.n_atoms, 3);

  for(int k = 0; k < 1; k++){
    /* Finding the optimal rotation 
     * In order to grant to find the right alignment
     * we do max 20 iterations!
     */
    for(int i = 0; i < m1.n_atoms; i++){
      A->data[i][0] = m1.atoms[i].coord.x;
      A->data[i][1] = m1.atoms[i].coord.y;
      A->data[i][2] = m1.atoms[i].coord.z;

      B->data[i][0] = m2.atoms[i].coord.x;
      B->data[i][1] = m2.atoms[i].coord.y;
      B->data[i][2] = m2.atoms[i].coord.z;
    }

    initMatrix(&R);
    initDVector(&t);
    rmsd = ICP(A, B, 1e-6, &R, &t, NULL);
    if(rmsd < best_rmsd){
      /* We apply the rotation to this rotation and we save as the best */
      ApplyTransformation(A, R, t, &A_best);
      best_rmsd = rmsd;
    }
    else{
      /*This method destroy the original rotation!*/
      RandomConformationRotation(&m1);
    }
    DelMatrix(&R);
    DelDVector(&t);
  }
  
  for(int i = 0; i < m1.n_atoms; i++){
    m1.atoms[i].coord.x = A_best->data[i][0];
    m1.atoms[i].coord.y = A_best->data[i][1];
    m1.atoms[i].coord.z = A_best->data[i][2];
  }
  
  DelMatrix(&A_best);
  DelMatrix(&A);
  DelMatrix(&B);
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
