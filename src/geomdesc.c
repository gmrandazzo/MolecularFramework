/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <math.h>
#include <scientific.h>

#include "atomanalysis.h"
#include "geomdesc.h"
#include "molecule.h"
#include "misc.h"
#include "periodic_table.h"
#include "mol3Dfields.h"


double GetDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return sqrt(square(x1-x2) + square(y1-y2) + square(z1-z2));
}

double GetDistance_(POINT p1, POINT p2)
{
  return sqrt(square(p1.x - p2.x) + square(p1.y - p2.y) + square(p1.z -p2.z));
}

/**
 * Calculates the angle (in radians) between two vectors pointing outward from one center
 *
 * p1 first point
 * p2 second point
 * c center point
 */
double GetAngle(POINT p1, POINT p2, POINT c)
{
  double p1c, p2c, p1p2, n, d, r;
  p1c = GetDistance_(p1, c);
  p2c = GetDistance_(p2, c);
  p1p2 = GetDistance_(p1, p2);
  n = (square(p1c) + square(p2c) - square(p1p2));
  d = (2*p1c*p2c);

  if(FLOAT_EQ(p1c, 0.f, 1e-5) || FLOAT_EQ(p1c, 0.f, 1e-5)){
    /*Atoms in the same plane*/
    return _pi_;
  }
  else{
    r = n/d;
    if(r > 1.f)
      return acos(1);
    else if (r < -1.f)
      return acos(-1);
    else
      return acos(n/d);
  }
}

double Rad2Grad(double rangle){
  return (rangle*360.)/ (2*_pi_);
}

double DihedralAngle(POINT p1, POINT p2, POINT p3, POINT p4)
{
  POINT q1, q2, q3, q12, q23;
  /*
  printf("Calculating dihedral angle for:\n");
  printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", p1.x, p1.y ,p1.z,
                                                                          p2.x, p2.y ,p2.z,
                                                                          p3.x, p3.y ,p3.z,
                                                                          p4.x, p4.y ,p4.z);
  */
  q1.x = (p2.x-p1.x);
  q1.y = (p2.y-p1.y);
  q1.z = (p2.z-p1.z);

  q2.x = (p3.x-p2.x);
  q2.y = (p3.y-p2.y);
  q2.z = (p3.z-p2.z);

  q3.x = (p4.x-p3.x);
  q3.y = (p4.y-p3.y);
  q3.z = (p4.z-p3.z);

  /*cross products*/
  q12.x = q1.y*q2.z - q1.z*q2.y;
  q12.y = q1.z*q2.x - q1.x*q2.z;
  q12.z = q1.x*q2.y - q1.y*q2.x;

  q23.x = q2.y*q3.z - q2.z*q3.y;
  q23.y = q2.z*q3.x - q2.x*q3.z;
  q23.z = q2.x*q3.y - q2.y*q3.x;

  double q12norm =  sqrt(q12.x*q12.x + q12.y*q12.y + q12.z*q12.z);
  double q23norm =  sqrt(q23.x*q23.x + q23.y*q23.y + q23.z*q23.z);

  if(FLOAT_EQ(q12norm, 0.f, 1e-5) || FLOAT_EQ(q23norm, 0.f, 1e-5)){
    /*The angle is 180 or the angle is not possible*/
    return _pi_;
  }
  else{
    q12.x /= q12norm;
    q12.y /= q12norm;
    q12.z /= q12norm;

    q23.x /= q23norm;
    q23.y /= q23norm;
    q23.z /= q23norm;

    double a = (q12.x*q3.x)+(q12.y*q3.y)+(q12.z*q3.z);

    int porm;
    if(a < 0)
      porm = -1;
    else if(FLOAT_EQ(a, 0.f, 1e-8))
      porm = 0;
    else
      porm = 1;

    double num = ((q12.x*q23.x)+(q12.y*q23.y)+(q12.z*q23.z));
    double den = sqrt(((q12.x*q12.x)+(q12.y*q12.y)+(q12.z*q12.z)) * ((q23.x*q23.x)+(q23.y*q23.y)+(q23.z*q23.z)));
    double rad = num/den;
    rad = (double)round(rad*100000.0 )/100000.0;
    if(rad > 1.f)
      rad = acos(1.f);
    else if(rad < -1.f)
      rad = acos(-1.f);
    else
      rad = acos(rad);

    if(porm > 0 || porm < 0)
      rad *= porm;
    return rad;
  }
}

/* Calculate the box size which contain the entire molecule.
 * This function is used by mol3DFields and shapedesc
 */
void GetMolBox(MOLECULE molecule, double size, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax)
{
  int l;
  double cx, cy, cz, hsize = size/2.f;
  cx = cy = cz = 0.f;
  for(l = 0; l < molecule.n_atoms; l++){
    cx += molecule.atoms[l].coord.x;
    cy += molecule.atoms[l].coord.y;
    cz += molecule.atoms[l].coord.z;
  }

  cx /= (double)molecule.n_atoms;
  cy /= (double)molecule.n_atoms;
  cz /= (double)molecule.n_atoms;

  (*xmin) = cx - hsize;
  (*xmax) = cx + hsize;
  (*ymin) = cy - hsize;
  (*ymax) = cy + hsize;
  (*zmin) = cz - hsize;
  (*zmax) = cz + hsize;
}

void GetPlanarity(MOLECULE molecule, double *planarity)
{
  size_t i;
  matrix *coord;

  if(molecule.n_atoms > 1){
    NewMatrix(&coord, molecule.n_atoms, 3);
    for(i = 0; i < molecule.n_atoms; i++){
      setMatrixValue(coord, i, 0, molecule.atoms[i].coord.x);
      setMatrixValue(coord, i, 1, molecule.atoms[i].coord.y);
      setMatrixValue(coord, i, 2, molecule.atoms[i].coord.z);
    }

    (*planarity) = CalcPlanarity(coord);

    DelMatrix(&coord);
  }
  else{
    (*planarity) = 0.f;
  }
}

double CalcPlanarity(matrix *coord)
{
  double planarity;
  matrix *residuals;
  dvector *colsdev;
  PCAMODEL *m;

  NewPCAModel(&m);
  PCA(coord, 0, 2, m, NULL);

  /* >>> compute planarity */
  initMatrix(&residuals);

  GetResidualMatrix(coord, m, 2, &residuals);

  /* Calculate STDEV for residuals and then the results is the descriptor */
  initDVector(&colsdev);
  MatrixColSDEV(residuals, &colsdev);

  planarity = sqrt(square(getDVectorValue(colsdev, 0))+square(getDVectorValue(colsdev, 1))+square(getDVectorValue(colsdev, 2)));

  DelMatrix(&residuals);
  DelDVector(&colsdev);
  DelPCAModel(&m);

  if(_isnan_(planarity)){
    return 0.f;
  }
  else{
    return planarity;
  }
}


void GetMoleculaLenght(MOLECULE molecule, double* lenght)
{
  size_t i, natom;
  matrix *coord;
  PCAMODEL *m;
  double t1min, t1max, t2min, t2max, l1, l2;

  natom = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    if(strcmp(molecule.atoms[i].asymbl, "H") != 0){ /* no hydrogen atoms */
      natom++;
    }
  }

  NewMatrix(&coord, natom, 3);
  natom = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    if(strcmp(molecule.atoms[i].asymbl, "H") != 0){
      setMatrixValue(coord, natom, 0, molecule.atoms[i].coord.x);
      setMatrixValue(coord, natom, 1, molecule.atoms[i].coord.y);
      setMatrixValue(coord, natom, 2, molecule.atoms[i].coord.z);
      natom++;
    }
  }

  NewPCAModel(&m);
  PCA(coord, 0, 2, m, NULL);

  /* >>> compute lenght which is the minimum point and the maximum point of the scores */
  MatrixColumnMinMax(m->scores, 0, &t1min, &t1max);
  MatrixColumnMinMax(m->scores, 1, &t2min, &t2max);

  l1 = sqrt(square(t1max-t1min));
  l2 = sqrt(square(t2max-t2min));

  if(l1 < l2)
    (*lenght) = l2;
  else
    (*lenght) = l1;

  DelPCAModel(&m);
  DelMatrix(&coord);
}


/*FP array position*/
#define MAX_3DFP_ATOMS 26
#define TYPEDIM 32

struct atom_types{
  char name[TYPEDIM];
  int pos;
};

static struct atom_types fp_atom_tab[MAX_3DFP_ATOMS] = {
  {"H", 0},
  {"C.1", 1},
  {"C.2", 2},
  {"C.3", 3},
  {"C.ar", 4},
  {"C.cat", 5},
  {"N.1", 6},
  {"N.2", 7},
  {"N.3", 8},
  {"N.4", 9},
  {"N.am", 10},
  {"N.ar", 11},
  {"N.pl3", 12},
  {"O.2", 13},
  {"O.3", 14},
  {"O.co2", 15},
  {"P.3", 16},
  {"S.2", 17},
  {"S.3", 18},
  {"S.O", 19},
  {"S.O2", 20},
  {"F", 21},
  {"Cl", 22},
  {"Br", 23},
  {"I", 24},
  {"Se", 25}
  /*
  {"Li", 5},
  {"Na", 5},
  {"K", 5},
  {"Be", 5},
  {"Mg", 5},
  {"Ca", 5},
  {"B", 5},
  {"Al", 5},
  {"Si", 5},
  {"Fe", 5},
  {"Co", 5},
  {"Zn", 5},
  {"Cu", 5},
  {"Mn", 5},
  */
};

int getAtomTypeFPPos(const char *type)
{
  int i, pos = -1;
  for(i = 0; i < MAX_3DFP_ATOMS; i++){
    if(strstr(type, fp_atom_tab[i].name) != NULL){
      pos = fp_atom_tab[i].pos;
      break;
    }
    else{
      continue;
    }
  }
  return pos;
}

/*
 * Convert 3D Fingerprint to a 2D matrix representation.
 * This representation is useful for convolutional neural networks
 */
void FP2Matrix(dvector *fp, matrix **m)
{
  size_t i, j, c;
  size_t nsteps = (fp->size/MAX_3DFP_ATOMS);
  ResizeMatrix(m, MAX_3DFP_ATOMS, nsteps);
  c = 0;
  for(i = 0; i < nsteps; i++){
    for(j = 0; j < MAX_3DFP_ATOMS; j++){
      (*m)->data[j][i] = fp->data[c];
      c++;
    }
  }
}

/* This methodology is ok and was already studied in 2D
 * see: https://www.researchgate.net/publication/263355029_Cytochrome_P450_site_of_metabolism_prediction_from_2D_topological_fingerprints_using_GPU_accelerated_probabilistic_classifiers/figures?lo=1
 */
void Get3DAtomSumDistanceFingerprint(MOLECULE molecule,
                                     int atom_id,
                                     double sph_step,
                                     double max_radius,
                                     dvector **fp)
{
  POINT from;
  from =  molecule.atoms[atom_id].coord;
  Get3DPosSumDistanceFingerprint(molecule, from, sph_step, max_radius, fp);
}

/* This methodology is ok and was already studied in 2D
 * see: https://www.researchgate.net/publication/263355029_Cytochrome_P450_site_of_metabolism_prediction_from_2D_topological_fingerprints_using_GPU_accelerated_probabilistic_classifiers/figures?lo=1
 */
void Get3DTwoAtomSumDistanceFingerprint(MOLECULE molecule,
                                     int atom_id1,
                                     int atom_id2,
                                     double sph_step,
                                     double max_radius,
                                     dvector **fp)
{
  POINT p;
  p =  molecule.atoms[atom_id1].coord;
  p.x += molecule.atoms[atom_id2].coord.x;
  p.y += molecule.atoms[atom_id2].coord.y;
  p.z += molecule.atoms[atom_id2].coord.z;

  p.x /= 2.;
  p.y /= 2.;
  p.z /= 2.;

  Get3DPosSumDistanceFingerprint(molecule, p, sph_step, max_radius, fp);
}

void Get3DPosSumDistanceFingerprint(MOLECULE molecule, POINT from, double sph_step, double max_radius, dvector **fp)
{
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT to;
  /*from.x = molecule.atoms[atom_id].coord.x;
  from.y = molecule.atoms[atom_id].coord.y;
  from.z = molecule.atoms[atom_id].coord.z;
  */
  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        (*fp)->data[pos+mult] += dst;
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}

/*
 * Potential sum fingerprint. This potential is generated
 * according to an atom type property "AtomsProperty *w_atypes"
 * These w_atypes are learned using an optimisation algorithm.
 */
void Get3DAtomSumPotentialFingerprint(MOLECULE molecule,
                                      int atom_id,
                                      double sph_step,
                                      double max_radius,
                                      AtomsProperty *w_atypes,
                                      dvector **fp)
{
  POINT from;
  from =  molecule.atoms[atom_id].coord;
  Get3DPosSumPotentialFingerprint(molecule, from, sph_step, max_radius, w_atypes, fp);
}

void Get3DTwoAtomSumPotentialFingerprint(MOLECULE molecule,
                                         int atom_id1,
                                         int atom_id2,
                                         double sph_step,
                                         double max_radius,
                                         AtomsProperty *w_atypes,
                                         dvector **fp)
{
  POINT p;
  p =  molecule.atoms[atom_id1].coord;
  p.x += molecule.atoms[atom_id2].coord.x;
  p.y += molecule.atoms[atom_id2].coord.y;
  p.z += molecule.atoms[atom_id2].coord.z;

  p.x /= 2.;
  p.y /= 2.;
  p.z /= 2.;

  Get3DPosSumPotentialFingerprint(molecule, p, sph_step, max_radius, w_atypes, fp);
}

void Get3DPosSumPotentialFingerprint(MOLECULE molecule,
                                     POINT from,
                                     double sph_step,
                                     double max_radius,
                                     AtomsProperty *w_atypes,
                                     dvector **fp)
{
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT to;
  /*from.x = molecule.atoms[atom_id].coord.x;
  from.y = molecule.atoms[atom_id].coord.y;
  from.z = molecule.atoms[atom_id].coord.z;
  */
  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        double weight = getGenericAProperty(molecule.atoms[i].type, w_atypes);
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        (*fp)->data[pos+mult] += weight/square(dst);
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}

/*
 * Scan over all atoms and find all the atoms attached to atom_id with in two
 * bond lengths a-p-b
 */
void GetAngleSumFingerprint(MOLECULE molecule, int atom_id, dvector **fp)
{
  size_t i, pos;
  double alpha;
  uivector *aid;

  /*
   int disc_angle = 5;
   from 100 degrees to 180 degrees
   * 180 = carbon dioxide
   * 120 = BF3
   * 109.5 = CH4
   * 104 = H2O
   *
   * 100/106 106/110 110/120  120/140 140/180
   */

  NewDVector(fp, MAX_3DFP_ATOMS);
  for(i = 0; i < molecule.n_atoms; i++){
    if(i != atom_id){
      initUIVector(&aid);
      Get3DAngle(molecule, atom_id, i, &alpha, &aid);
      if(FLOAT_EQ(alpha, -9999.0, 1e-1)){
        DelUIVector(&aid);
        continue;
      }
      else{
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        (*fp)->data[pos] += fabs(alpha);
        /*printf("angle %zu %s  %zu %s  %zu %s -> %f\n", aid->data[0],
                                                       molecule.atoms[aid->data[0]].type,
                                                       aid->data[1],
                                                       molecule.atoms[aid->data[1]].type,
                                                       aid->data[2],
                                                       molecule.atoms[aid->data[2]].type,
                                                       alpha);*/

      }
      DelUIVector(&aid);
    }
    else{
      continue;
    }
  }
}

/*
 * Scan over all atoms and find all the atoms attached to atom_id with in three
 * bond lengths a-p1-p2-b
 */
void GetTorsionSumFingerprint(MOLECULE molecule, int atom_id, dvector **fp)
{
  size_t i, pos;
  double delta;
  NewDVector(fp, MAX_3DFP_ATOMS);
  for(i = 0; i < molecule.n_atoms; i++){
    if(i != atom_id){
     Get3DDihedralAngle(molecule, atom_id, i, &delta, NULL);
     if(FLOAT_EQ(delta, -9999.0, 1e-1)){
       continue;
      }
      else{
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        (*fp)->data[pos] += fabs(delta);
      }
    }
    else{
      continue;
    }
  }
}


/*Here start the code with a bit of fantasy....*/

void Get3DBaricentreSumAngleAtomFingerprint(MOLECULE molecule,
                                            int atom_id,
                                            double sph_step,
                                            double max_radius,
                                            dvector **fp)
{
  POINT from;
  from =  molecule.atoms[atom_id].coord;
  Get3DBaricentreSumAngleFingerprint(molecule, from, sph_step, max_radius, fp);
}


void Get3DBaricentreSumAngleTwoAtomFingerprint(MOLECULE molecule,
                                               int atom_id1,
                                               int atom_id2,
                                               double sph_step,
                                               double max_radius,
                                               dvector **fp)
{
  POINT p;
  p =  molecule.atoms[atom_id1].coord;
  p.x += molecule.atoms[atom_id2].coord.x;
  p.y += molecule.atoms[atom_id2].coord.y;
  p.z += molecule.atoms[atom_id2].coord.z;

  p.x /= 2.;
  p.y /= 2.;
  p.z /= 2.;

  Get3DBaricentreSumAngleFingerprint(molecule, p, sph_step, max_radius, fp);
}


void Get3DBaricentreSumAngleFingerprint(MOLECULE molecule,
                                        POINT from,
                                        double sph_step,
                                        double max_radius,
                                        dvector **fp)
{
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT to, baricentre;
  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*printf("baricenter -> x: %f y: %f z: %f\n", baricentre.x, baricentre.y, baricentre.z);*/
  /*from.x = molecule.atoms[atom_id].coord.x;
  from.y = molecule.atoms[atom_id].coord.y;
  from.z = molecule.atoms[atom_id].coord.z;
  */
  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        /* Calculate the angle */
        (*fp)->data[pos+mult] += Rad2Grad(GetAngle(from, to, baricentre));
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}


void Get3DAtomBaricenterAtomMagnitudeFingerprint(MOLECULE molecule,
                                                 int atom_id,
                                                 double sph_step,
                                                 double max_radius,
                                                 dvector **fp)
{
  POINT from;
  from =  molecule.atoms[atom_id].coord;
  Get3DAtomBaricenterMagnitudeFingerprint(molecule, from, sph_step, max_radius, fp);
}


void Get3DTwoAtomBaricenterAtomMagnitudeFingerprint(MOLECULE molecule,
                                                    int atom_id1,
                                                    int atom_id2,
                                                    double sph_step,
                                                    double max_radius,
                                                    dvector **fp)
{
  POINT p;
  p =  molecule.atoms[atom_id1].coord;
  p.x += molecule.atoms[atom_id2].coord.x;
  p.y += molecule.atoms[atom_id2].coord.y;
  p.z += molecule.atoms[atom_id2].coord.z;

  p.x /= 2.;
  p.y /= 2.;
  p.z /= 2.;

  Get3DBaricentreSumAngleFingerprint(molecule, p, sph_step, max_radius, fp);
}


void Get3DAtomBaricenterMagnitudeFingerprint(MOLECULE molecule,
                                             POINT from,
                                             double sph_step,
                                             double max_radius,
                                             dvector **fp)
{
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT to, baricentre;
  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*printf("baricenter -> x: %f y: %f z: %f\n", baricentre.x, baricentre.y, baricentre.z);*/
  /*from.x = molecule.atoms[atom_id].coord.x;
  from.y = molecule.atoms[atom_id].coord.y;
  from.z = molecule.atoms[atom_id].coord.z;
  */
  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        pos = getAtomTypeFPPos(molecule.atoms[i].type);
        /* Calculate the angle
        (*fp)->data[pos+mult] += GetAngle(from, to, baricentre);
        */
        /* Calcualate the dipole moment between From-Baricenter-TO */
        double mux = from.x + baricentre.x + to.x;
        double muy = from.y + baricentre.y + to.y;
        double muz = from.z + baricentre.z + to.z;
        (*fp)->data[pos+mult] += sqrt(square(mux)+square(muy)+square(muz));
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}

void Get3DTwoAtomBaricentreSumDihedralAngleFingerprint(MOLECULE molecule,
                                                      int atom_id1,
                                                      int atom_id2,
                                                      double sph_step,
                                                      double max_radius,
                                                      dvector **fp)
{
  /* dihedral angle calculated between p1, p2, p3, p4
   * p1 and p2 are atom_id1 and atom_id2 arranged to be the nearest neighbour to
   * the atom found in the sphere.
   * p3 is always the molecule baricentre
   * p4 is the atom found in sphere
   */
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT a1, a2, from, to, baricentre;
  a1 = molecule.atoms[atom_id1].coord;
  a2 = molecule.atoms[atom_id2].coord;


  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*From is utilised to build the sphere fingerprint
   * It is a medium point to check wich atoms are in this distance
   */
  from =  molecule.atoms[atom_id1].coord;
  from.x += molecule.atoms[atom_id2].coord.x;
  from.y += molecule.atoms[atom_id2].coord.y;
  from.z += molecule.atoms[atom_id2].coord.z;

  from.x /= 2.;
  from.y /= 2.;
  from.z /= 2.;

  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        pos = getAtomTypeFPPos(molecule.atoms[i].type);

        double a1_x_to = GetDistance_(a1, to);
        double a2_x_to = GetDistance_(a2, to);
        if(FLOAT_EQ(a1_x_to, 0.f, 1e-3) || FLOAT_EQ(a2_x_to, 0.f, 1e-3)){
          continue;
        }
        else{
          /*
          printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", a1.x, a1.y, a1.z,
                                                              a2.x, a2.y, a2.z,
                                                              baricentre.x, baricentre.y, baricentre.z,
                                                              to.x, to.y, to.z,
                                                              DihedralAngle(a1, a2, baricentre, to),
                                                              DihedralAngle(a2, a1, baricentre, to));
          */
          /*
          (*fp)->data[pos+mult] += fabs(Rad2Grad(DihedralAngle(a1, a2, baricentre, to)));
          (*fp)->data[pos+mult] += fabs(Rad2Grad(DihedralAngle(a2, a1, baricentre, to)));
          */
          (*fp)->data[pos+mult] += fabs(DihedralAngle(a1, a2, baricentre, to));
          (*fp)->data[pos+mult] += fabs(DihedralAngle(a2, a1, baricentre, to));
        }
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}


void Get3DTwoAtomDistanceDihedralAngleCombinationFingerprint(MOLECULE molecule,
                                                             int atom_id1,
                                                             int atom_id2,
                                                             double sph_step,
                                                             double max_radius,
                                                             dvector **fp)
{
  /*
   * the distance is calculate from the medium point.
   *
   * dihedral angle calculated between p1, p2, p3, p4
   * p1 and p2 are atom_id1 and atom_id2 arranged to be the nearest neighbour to
   * the atom found in the sphere.
   * p3 is always the molecule baricentre
   * p4 is the atom found in sphere
   */
  size_t s, i, pos, mult;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f;
  POINT a1, a2, from, to, baricentre;
  a1 = molecule.atoms[atom_id1].coord;
  a2 = molecule.atoms[atom_id2].coord;


  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*From is utilised to build the sphere fingerprint
   * It is a medium point to check wich atoms are in this distance
   */
  from =  molecule.atoms[atom_id1].coord;
  from.x += molecule.atoms[atom_id2].coord.x;
  from.y += molecule.atoms[atom_id2].coord.y;
  from.z += molecule.atoms[atom_id2].coord.z;

  from.x /= 2.;
  from.y /= 2.;
  from.z /= 2.;

  NewDVector(fp, nsteps*MAX_3DFP_ATOMS);
  mult = 0;
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        pos = getAtomTypeFPPos(molecule.atoms[i].type);

        double a1_x_to = GetDistance_(a1, to);
        double a2_x_to = GetDistance_(a2, to);
        if(FLOAT_EQ(a1_x_to, 0.f, 1e-3) || FLOAT_EQ(a2_x_to, 0.f, 1e-3)){
          continue;
        }
        else{
          /*
          printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", a1.x, a1.y, a1.z,
                                                              a2.x, a2.y, a2.z,
                                                              baricentre.x, baricentre.y, baricentre.z,
                                                              to.x, to.y, to.z,
                                                              DihedralAngle(a1, a2, baricentre, to),
                                                              DihedralAngle(a2, a1, baricentre, to));
          */

          /*influence of atom "to" into atom "a1"*/
          (*fp)->data[pos+mult] += dst*fabs(DihedralAngle(a1, a2, baricentre, to));
          /*influence of atom "to" into atom "a2"*/
          (*fp)->data[pos+mult] += dst*fabs(DihedralAngle(a2, a1, baricentre, to));
        }
      }
      else{
        continue;
      }
    }
    x+=sph_step;
    mult+=MAX_3DFP_ATOMS;
  }
}

void Get3DEPotDihedralAngleWeightedAtomFingerprint(MOLECULE molecule,
                                                   int atom_id,
                                                   double sph_step,
                                                   double max_radius,
                                                   dvector **fp)
{
  /* For each radius see which atom is in radius.
   * for each atom in radius calculate the elctrostatic potential
   * and rescale this potential for the dihedral angle formed with the
   * medium point and baricenter
   * dihedral angle calculated between p1, p2, p3, p4
   * p1 and p2 are atom_id and is connection arranged to be
   * the nearest neighbour of the baricentre
   * p3 is always the molecule baricentre
   * p4 is the atom found in sphere
   */
  size_t s, i, j;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f, el_pot;
  double dst = 0.f;
  POINT from, p, to, baricentre;

  /* Atom a is bonded with */
  uivector *a_conn;
  initUIVector(&a_conn);
  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].origin_atom_id == atom_id){
      UIVectorAppend(&a_conn, molecule.bonds[i].target_atom_id);
    }
    else if(molecule.bonds[i].target_atom_id == atom_id){
      UIVectorAppend(&a_conn, molecule.bonds[i].origin_atom_id);
    }
    else{
      continue;
    }
  }

  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*From is utilised to build the sphere fingerprint
   * It is a medium point to check wich atoms are in this distance
   */
  from =  molecule.atoms[atom_id].coord;

  NewDVector(fp, nsteps);
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        /* calculate the Electrostatic potential value at this distance
         * using this atom and find the position were you should count the value.
         */
        dst *= 1E-10; /*Convert the distance from angstrom to meter*/
        /* calculate the potential in Joules */
        el_pot = (epsilon0_C_J_m*au2C(molecule.atoms[i].charge))/dst;
        /* Conversion from J to kcal/mol dividing the value for
        * 4184 J / 6.022E23 mol (avogadro constant)
        */
        el_pot /= 6.947857854533377e-21;
        // printf("el_pot: %e %e\n", dst, el_pot);

        double res = 0.f;
        for(j = 0; j < a_conn->size; j++){
          p = molecule.atoms[a_conn->data[i]].coord;
          double dh = DihedralAngle(from, p, baricentre, to);
          double k;
          if(FLOAT_EQ(dh, -9999., 1e-3)){
            /*
             * Angle not possible... this means no dihedral effect
             */
             k = 0.f;
          }
          else{
            double delta = fabs(Rad2Grad(dh));
            /*Karplus equation weight*/
            k = cos(2*delta)+cos(delta);
          }

          // printf("KARPLUS %f %f\n", k);
          res += el_pot*k;
        }

        (*fp)->data[s] += res/(double)a_conn->size;
      }
      else{
        continue;
      }
    }
    x+=sph_step;
  }
}

void Get3DEPotDihedralAngleWeightedTwoAtomFingerprint(MOLECULE molecule,
                                                      int atom_id1,
                                                      int atom_id2,
                                                      double sph_step,
                                                      double max_radius,
                                                      dvector **fp)
{
  /* For each radius see which atom is in radius.
   * for each atom in radius calculate the elctrostatic potential
   * and rescale this potential for the dihedral angle formed with the
   * medium point and baricenter
   * dihedral angle calculated between p1, p2, p3, p4
   * p1 and p2 are atom_id1 and atom_id2 arranged to be the nearest neighbour to
   * the atom found in the sphere.
   * p3 is always the molecule baricentre
   * p4 is the atom found in sphere
   */
  size_t s, i;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f, el_pot;
  double dst = 0.f;
  POINT a1, a2, from, to, baricentre;
  a1 = molecule.atoms[atom_id1].coord;
  a2 = molecule.atoms[atom_id2].coord;


  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*From is utilised to build the sphere fingerprint
   * It is a medium point to check wich atoms are in this distance
   */
  from =  molecule.atoms[atom_id1].coord;
  from.x += molecule.atoms[atom_id2].coord.x;
  from.y += molecule.atoms[atom_id2].coord.y;
  from.z += molecule.atoms[atom_id2].coord.z;

  from.x /= 2.;
  from.y /= 2.;
  from.z /= 2.;

  NewDVector(fp, nsteps);
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        /* calculate the Electrostatic potential value at this distance
         * using this atom and find the position were you should count the value.
         */
        dst *= 1E-10; /*Convert the distance from angstrom to meter*/
        /* calculate the potential in Joules */
        el_pot = (epsilon0_C_J_m*au2C(molecule.atoms[i].charge))/dst;
        /* Conversion from J to kcal/mol dividing the value for
        * 4184 J / 6.022E23 mol (avogadro constant)
        */
        el_pot /= 6.947857854533377e-21;
        // printf("el_pot: %e %e\n", dst, el_pot);

        double a1_x_to = GetDistance_(a1, to);
        double a2_x_to = GetDistance_(a2, to);
        if(FLOAT_EQ(a1_x_to, 0.f, 1e-3) || FLOAT_EQ(a2_x_to, 0.f, 1e-3)){
          continue;
        }
        else{
          /*
          printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", a1.x, a1.y, a1.z,
                                                              a2.x, a2.y, a2.z,
                                                              baricentre.x, baricentre.y, baricentre.z,
                                                              to.x, to.y, to.z,
                                                              DihedralAngle(a1, a2, baricentre, to),
                                                              DihedralAngle(a2, a1, baricentre, to));
          */

          double dh1 = DihedralAngle(a1, a2, baricentre, to);
          double dh2 = DihedralAngle(a2, a1, baricentre, to);
          double k1, k2;
          if(FLOAT_EQ(dh1, -9999., 1e-3)){
            /*
             * Angle not possible... this means no dihedral effect
             */
             k1 = 0.f;
          }
          else{
            double delta1 = fabs(Rad2Grad(dh1));
            /*Karplus equation weight*/
            k1 = cos(2*delta1)+cos(delta1);
          }

          if(FLOAT_EQ(dh2, -9999., 1e-3)){
            /*
             * Angle not possible... this means no dihedral effect
             */
            k2 = 0.f;
          }
          else{
            double delta2 = fabs(Rad2Grad(dh2));
            k2 = cos(2*delta2)+cos(delta2);
          }

          // printf("KARPLUS %f %f\n", k1, k2);
          (*fp)->data[s] += el_pot*k1;
          (*fp)->data[s] += el_pot*k2;
        }
      }
      else{
        continue;
      }
    }
    x+=sph_step;
  }
}


void Get3DEPotAngleWeightedAtomFingerprint(MOLECULE molecule,
                                           int atom_id,
                                           double sph_step,
                                           double max_radius,
                                           dvector **fp)
{
  POINT from;
  from =  molecule.atoms[atom_id].coord;
  Get3DEPotAngleWeightedFingerprint(molecule, from, sph_step, max_radius, fp);
}

void Get3DEPotAngleWeightedTwoAtomFingerprint(MOLECULE molecule,
                                              int atom_id1,
                                              int atom_id2,
                                              double sph_step,
                                              double max_radius,
                                              dvector **fp)
{
  POINT from;

  from =  molecule.atoms[atom_id1].coord;
  from.x += molecule.atoms[atom_id2].coord.x;
  from.y += molecule.atoms[atom_id2].coord.y;
  from.z += molecule.atoms[atom_id2].coord.z;

  from.x /= 2.;
  from.y /= 2.;
  from.z /= 2.;
  Get3DEPotAngleWeightedFingerprint(molecule, from, sph_step, max_radius, fp);
}

void Get3DEPotAngleWeightedFingerprint(MOLECULE molecule,
                                       POINT from,
                                       double sph_step,
                                       double max_radius,
                                       dvector **fp)
{
  size_t s, i;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst = 0.f, el_pot;
  POINT to, baricentre;

  baricentre.x = baricentre.y = baricentre.z = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    baricentre.x += molecule.atoms[i].coord.x;
    baricentre.y += molecule.atoms[i].coord.y;
    baricentre.z += molecule.atoms[i].coord.z;
  }

  baricentre.x /= (double)molecule.n_atoms;
  baricentre.y /= (double)molecule.n_atoms;
  baricentre.z /= (double)molecule.n_atoms;

  /*printf("baricenter -> x: %f y: %f z: %f\n", baricentre.x, baricentre.y, baricentre.z);*/
  /*from.x = molecule.atoms[atom_id].coord.x;
  from.y = molecule.atoms[atom_id].coord.y;
  from.z = molecule.atoms[atom_id].coord.z;
  */
  NewDVector(fp, nsteps);
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      /*to.x = molecule.atoms[i].coord.x;
      to.y = molecule.atoms[i].coord.y;
      to.z = molecule.atoms[i].coord.z;*/
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if(dst > x && dst < x+sph_step){
        /* Calculate the angle */
        dst *= 1E-10; /*Convert the distance from angstrom to meter*/
        /* calculate the potential in Joules */
        el_pot = (epsilon0_C_J_m*au2C(molecule.atoms[i].charge))/dst;
        /* Conversion from J to kcal/mol dividing the value for
        * 4184 J / 6.022E23 mol (avogadro constant)
        */
        el_pot /= 6.947857854533377e-21;
        //printf("el_pot: %e %e\n", dst, el_pot);
        double angle = GetAngle(from, to, baricentre);
        if(FLOAT_EQ(angle, -9999., 1e-3)){
          /*No angle effect possible*/
          continue;
        }
        else{
          (*fp)->data[s] += el_pot*Rad2Grad(angle);
        }
      }
      else{
        continue;
      }
    }
    x+=sph_step;
  }
}

/* Here the code OK */
void Get3DAtomEpotFingerprint(MOLECULE molecule,
                              int atom_id,
                              double sph_step,
                              double max_radius,
                              dvector **fp)
{
  POINT from;
  from = molecule.atoms[atom_id].coord;
  Get3DEpotFingerprint(molecule, from, sph_step, max_radius, fp);
}


void Get3DTwoAtomEpotFingerprint(MOLECULE molecule,
                                 int atom_id1,
                                 int atom_id2,
                                 double sph_step,
                                 double max_radius,
                                 dvector **fp)
{
  POINT from;

  from = molecule.atoms[atom_id1].coord;
  from.x += molecule.atoms[atom_id2].coord.x;
  from.y += molecule.atoms[atom_id2].coord.y;
  from.z += molecule.atoms[atom_id2].coord.z;
  from.x /= 2.f;
  from.y /= 2.f;
  from.z /= 2.f;

  Get3DEpotFingerprint(molecule, from, sph_step, max_radius, fp);
}

/* This methodology is ok and was already studied in 2D
 * see: https://www.researchgate.net/publication/263355029_Cytochrome_P450_site_of_metabolism_prediction_from_2D_topological_fingerprints_using_GPU_accelerated_probabilistic_classifiers/figures?lo=1
 */
void Get3DEpotFingerprint(MOLECULE molecule,
                          POINT from,
                          double sph_step,
                          double max_radius,
                          dvector **fp)
{
  size_t s, i;
  size_t nsteps = (size_t)ceil(max_radius/sph_step);
  double x = 0.f;
  double dst;

  double el_pot = 0.f;
  POINT to;

  NewDVector(fp, nsteps);
  for(s = 0; s < nsteps; s++){
    for(i = 0; i < molecule.n_atoms; i++){
      to = molecule.atoms[i].coord;
      dst = GetDistance_(from, to);
      if((dst > x && dst < x+sph_step)){
        /* calculate the Electrostatic potential value at this distance
         * using this atom and find the position were you should count the value.
         */
        dst *= 1E-10; /*Convert the distance from angstrom to meter*/
        /* calculate the potential in Joules */
        el_pot = (epsilon0_C_J_m*au2C(molecule.atoms[i].charge))/dst;
        /* Conversion from J to kcal/mol dividing the value for
        * 4184 J / 6.022E23 mol (avogadro constant)
        */
        el_pot /= 6.947857854533377e-21;
        //printf("el_pot: %e %e\n", dst, el_pot);
        (*fp)->data[s] += el_pot;
      }
    }
    x+=sph_step;
  }
}

/*Function utilised by Get3DDihedralAngle */
int is_p1abp2(MOLECULE molecule, int ai, int aj)
{
  int k;
  for(k = 0; k < molecule.n_bonds; k++){
    if((ai == molecule.bonds[k].origin_atom_id &&
        aj  == molecule.bonds[k].target_atom_id) ||
      (ai == molecule.bonds[k].target_atom_id &&
       aj  == molecule.bonds[k].origin_atom_id)){
      // this is the couple p1-a-b-p2;
      return 1;
    }
    else{
      continue;
    }
  }
  return 0;
}

void Get3DDihedralAngle(MOLECULE molecule,
                        int atom_id1,
                        int atom_id2,
                        double *dh_angle,
                        uivector **dhids)
{
  size_t i, j, ai, aj;
  int compute_dihedral = 0;
  /* Find a and b wich satisfy the following condition
  * p1-a-b-p2
  */
  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].origin_atom_id == atom_id1){ // check if target is bonded with...
      for(j = 0; j < molecule.n_bonds; j++){
        if(i != j){
          if(molecule.bonds[j].origin_atom_id == atom_id2){
            if(is_p1abp2(molecule,
                         molecule.bonds[i].target_atom_id,
                         molecule.bonds[j].target_atom_id) == 1){
              ai = molecule.bonds[i].target_atom_id;
              aj = molecule.bonds[j].target_atom_id;
              compute_dihedral = 1;
              break;
            }
            else{
              continue;
            }
          }
          else if(molecule.bonds[j].target_atom_id == atom_id2){
            if(is_p1abp2(molecule,
                         molecule.bonds[i].target_atom_id,
                         molecule.bonds[j].origin_atom_id) == 1){
              ai = molecule.bonds[i].target_atom_id;
              aj = molecule.bonds[j].origin_atom_id;
              compute_dihedral = 1;
              break;
            }
            else{
              continue;
            }

          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }
    }
    else if(molecule.bonds[i].target_atom_id == atom_id1){ // then is origin
      for(j = 0; j < molecule.n_bonds; j++){
        if(i != j){
          if(molecule.bonds[j].origin_atom_id == atom_id2){
            if(is_p1abp2(molecule,
                         molecule.bonds[i].origin_atom_id,
                         molecule.bonds[j].target_atom_id) == 1){
              ai = molecule.bonds[i].origin_atom_id;
              aj = molecule.bonds[j].target_atom_id;
              compute_dihedral = 1;
              break;
            }
            else{
              continue;
            }
          }
          else if(molecule.bonds[j].target_atom_id == atom_id2){
            if(is_p1abp2(molecule,
                         molecule.bonds[i].origin_atom_id,
                         molecule.bonds[j].origin_atom_id) == 1){
              ai = molecule.bonds[i].origin_atom_id;
              aj = molecule.bonds[j].origin_atom_id;
              compute_dihedral = 1;
              break;
            }
            else{
              continue;
            }

          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }

    }
    else{
      continue;
    }
  }

  if(compute_dihedral == 1 &&
    atom_id1 != ai &&
    atom_id1 != aj &&
    atom_id2 != ai &&
    atom_id2 != aj &&
    atom_id1 != atom_id2 &&
    ai != aj){
    /*printf("%d %zu %zu %d\n", atom_id1+1, ai+1, aj+1, atom_id2+1);*/
    /*printf("%d %zu %zu %d\n", atom_id1, ai, aj, atom_id2);*/
    /* calculate the three vector q1 = atom_id1-ai, q2 = ai-aj, q3= aj-atom_id2*/
    (*dh_angle) = DihedralAngle(molecule.atoms[atom_id1].coord,
                                molecule.atoms[ai].coord,
                                molecule.atoms[aj].coord,
                                molecule.atoms[atom_id2].coord);
    if(_isnan_((*dh_angle))){
      (*dh_angle) = -9999.;
    }
    else{
      if(dhids != NULL){
        UIVectorResize(dhids, 4);
        (*dhids)->data[0] = atom_id1;
        (*dhids)->data[1] = ai;
        (*dhids)->data[2] = aj;
        (*dhids)->data[3] = atom_id2;
      }
    }
  }
  else{
    (*dh_angle) = -9999.;
  }
  return;
}


void Get3DAngle(MOLECULE molecule,
                int atom_id1,
                int atom_id2,
                double *angle,
                uivector **aid)
{
  size_t i, j, c;
  int compute_angle = 0;
  /* Find a the following condition
  * p1-a-p2
  */
  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].origin_atom_id == atom_id1){ // check if target is bonded with...
      for(j = 0; j < molecule.n_bonds; j++){
        if(i != j){
          if(molecule.bonds[j].origin_atom_id == atom_id2){
            if(molecule.bonds[i].target_atom_id == molecule.bonds[j].target_atom_id){
              // That is the angle which connect the two atoms...
              c = molecule.bonds[i].target_atom_id;
              compute_angle = 1;
              break;
            }
            else{
              continue;
            }
          }
          else if(molecule.bonds[j].target_atom_id == atom_id2){
            if(molecule.bonds[i].target_atom_id == molecule.bonds[j].origin_atom_id){
              // That is the angle which connect the two atoms...
              c = molecule.bonds[i].target_atom_id;
              compute_angle = 1;
              break;
            }
            else{
              continue;
            }
          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }
    }
    else if(molecule.bonds[i].target_atom_id == atom_id1){ // then is origin
      for(j = 0; j < molecule.n_bonds; j++){
        if(i != j){
          if(molecule.bonds[j].origin_atom_id == atom_id2){
            if(molecule.bonds[i].origin_atom_id == molecule.bonds[j].target_atom_id){
              // That is the angle which connect the two atoms...
              c = molecule.bonds[i].origin_atom_id;
              compute_angle = 1;
              break;
            }
            else{
              continue;
            }
          }
          else if(molecule.bonds[j].target_atom_id == atom_id2){
            if(molecule.bonds[i].origin_atom_id == molecule.bonds[j].origin_atom_id){
              // That is the angle which connect the two atoms...
              c = molecule.bonds[i].origin_atom_id;
              compute_angle = 1;
              break;
            }
            else{
              continue;
            }
          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }
    }
    else{
      continue;
    }
  }

  if(compute_angle == 1){
    //printf("%d %zu %zu %d\n", atom_id1+1, ai+1, aj+1, atom_id2+1);
    /* calculate the three vector q1 = atom_id1-ai, q2 = ai-aj, q3= aj-atom_id2*/
    (*angle) = GetAngle(molecule.atoms[atom_id1].coord,
                           molecule.atoms[atom_id2].coord,
                           molecule.atoms[c].coord);
    if(_isnan_((*angle)))
      (*angle) = -9999.;
    else{
      if(aid != NULL){
        UIVectorResize(aid, 3);
        (*aid)->data[0] = atom_id1;
        (*aid)->data[1] = c;
        (*aid)->data[2] = atom_id2;
      }
    }
  }
  else{
    (*angle) = -9999.;
  }
  return;
}
