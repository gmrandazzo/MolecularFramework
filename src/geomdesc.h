/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef GEOMDESC_H
#define GEOMDESC_H

#include <scientific.h>

#include "molecule.h"
#include "periodic_table.h"


/* Calculates the distance between two points */
double GetDistance(double x1, double y1, double z1, double x2, double y2, double z2);

/* Calculates the angle (in radians) between two vectors pointing outward from one center */
double GetDistance_(POINT p1, POINT p2);

/*Calculate angle between p1-c-p2 in radians */
double GetAngle(POINT p1, POINT p2, POINT c);

/*Calculate diedral angle in radians */
double DihedralAngle(POINT p1, POINT p2, POINT p3, POINT p4);

/*Convert angle from radians to grade */
double Rad2Grad(double rangle);

/* Calculate the box size which contain the entire molecule.
 * This function is used by mol3DFields and shapedesc
 */
void GetMolBox(MOLECULE molecule, double size, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax);

/*
 * Molecular Planarity
 */
void GetPlanarity(MOLECULE molecule, double *planarity);

/*
 * Calculate if some coordinates are planar or not
 */
double CalcPlanarity(matrix *coord);


/* Molecular Lenght
 *
 */
void GetMoleculaLenght(MOLECULE molecule, double *lenght);


/*
* Molecular diffusion using the stock einstein equation
* Algorithm described here: DOI: 10.1021/acs.molpharmaceut.7b01053
*/
void CalcMolecularDiffusion_StockesEinstein(MOLECULE molecule, double *diff);
void CalcMolecularDiffusion_WilkeChang(MOLECULE molecule, double *diff);

/*
 * Calculate Atomic Geometrical Fingerprint
 */

void FP2Matrix(dvector *fp, matrix **m);

void Get3DAtomSumDistanceFingerprint(MOLECULE molecule,
                                     int atom_id,
                                     double sph_step,
                                     double max_radius,
                                     dvector **fp);

void Get3DTwoAtomSumDistanceFingerprint(MOLECULE molecule,
                                        int atom_id1,
                                        int atom_id2,
                                        double sph_step,
                                        double max_radius,
                                        dvector **fp);

void Get3DPosSumDistanceFingerprint(MOLECULE molecule,
                                    POINT from,
                                    double sph_step,
                                    double max_radius,
                                    dvector **fp);

void Get3DAtomSumPotentialFingerprint(MOLECULE molecule,
                                      int atom_id,
                                      double sph_step,
                                      double max_radius,
                                      AtomsProperty *w_atypes,
                                      dvector **fp);

void Get3DTwoAtomSumPotentialFingerprint(MOLECULE molecule,
                                         int atom_id1,
                                         int atom_id2,
                                         double sph_step,
                                         double max_radius,
                                         AtomsProperty *w_atypes,
                                         dvector **fp);

void Get3DPosSumPotentialFingerprint(MOLECULE molecule,
                                    POINT from,
                                    double sph_step,
                                    double max_radius,
                                    AtomsProperty *w_atypes,
                                    dvector **fp);

void GetAngleSumFingerprint(MOLECULE molecule, int atom_id, dvector **fp);

void GetTorsionSumFingerprint(MOLECULE molecule, int atom_id, dvector **fp);

/*Here start the code with a bit of fantasy....*/
void Get3DBaricentreSumAngleAtomFingerprint(MOLECULE molecule,
                                         int atom_id,
                                         double sph_step,
                                         double max_radius,
                                         dvector **fp);

void Get3DBaricentreSumAngleTwoAtomFingerprint(MOLECULE molecule,
                                           int atom_id1,
                                           int atom_id2,
                                           double sph_step,
                                           double max_radius,
                                           dvector **fp);

void Get3DBaricentreSumAngleFingerprint(MOLECULE molecule,
                                        POINT from,
                                        double sph_step,
                                        double max_radius,
                                        dvector **fp);

void Get3DAtomBaricenterAtomMagnitudeFingerprint(MOLECULE molecule,
                                                 int atom_id,
                                                 double sph_step,
                                                 double max_radius,
                                                 dvector **fp);


void Get3DTwoAtomBaricenterAtomMagnitudeFingerprint(MOLECULE molecule,
                                                   int atom_id1,
                                                   int atom_id2,
                                                   double sph_step,
                                                   double max_radius,
                                                   dvector **fp);

void Get3DAtomBaricenterMagnitudeFingerprint(MOLECULE molecule,
                                            POINT from,
                                            double sph_step,
                                            double max_radius,
                                            dvector **fp);

void Get3DTwoAtomBaricentreSumDihedralAngleFingerprint(MOLECULE molecule,
                                                       int atom_id1,
                                                       int atom_id2,
                                                       double sph_step,
                                                       double max_radius,
                                                       dvector **fp);

void Get3DTwoAtomDistanceDihedralAngleCombinationFingerprint(MOLECULE molecule,
                                                             int atom_id1,
                                                             int atom_id2,
                                                             double sph_step,
                                                             double max_radius,
                                                             dvector **fp);

void Get3DEPotDihedralAngleWeightedAtomFingerprint(MOLECULE molecule,
                                                   int atom_id,
                                                   double sph_step,
                                                   double max_radius,
                                                   dvector **fp);

void Get3DEPotDihedralAngleWeightedTwoAtomFingerprint(MOLECULE molecule,
                                                      int atom_id1,
                                                      int atom_id2,
                                                      double sph_step,
                                                      double max_radius,
                                                      dvector **fp);

/*
 * Calculate the electrostatic potential effect on a specified atom
 * and weight this for the angle formed between atom, the atom found at distance d
 * and the baricenter of molecule
 */
void Get3DEPotAngleWeightedAtomFingerprint(MOLECULE molecule,
                                           int atom_id1,
                                           double sph_step,
                                           double max_radius,
                                           dvector **fp);

/*
* Calculate the electrostatic potential effect on two specified atoms
* and weight this for the angle formed between the medium point obtained by
* the x, y, z coordinates of the two atoms, the atom found at distance d
* and the baricenter of molecule
*/
void Get3DEPotAngleWeightedTwoAtomFingerprint(MOLECULE molecule,
                                             int atom_id1,
                                             int atom_id2,
                                             double sph_step,
                                             double max_radius,
                                             dvector **fp);

/*
 * Calculate the electrostatic potential effect on a specified point x, y, z
 * and weight this for the angle formed between the specified point x, y ,z
 * the atom found at distance d and the baricenter of molecule
 */
void Get3DEPotAngleWeightedFingerprint(MOLECULE molecule,
                                       POINT from,
                                       double sph_step,
                                       double max_radius,
                                       dvector **fp);


/*Here end the code with a bit of fantasy and start the code OK!!*/

/*
 * Calculate the electrostatic potential effect on a specified atom
 */
void Get3DAtomEpotFingerprint(MOLECULE molecule,
                              int atom_id,
                              double sph_step,
                              double max_radius,
                              dvector **fp);

/*
 * Calculate the electrostatic potential effect on two specified atoms
 */
void Get3DTwoAtomEpotFingerprint(MOLECULE molecule,
                                 int atom_id1,
                                 int atom_id2,
                                 double sph_step,
                                 double max_radius,
                                 dvector **fp);

/*
 * Calculate the electrostatic potential effect on a specified point x, y, z
 */
void Get3DEpotFingerprint(MOLECULE molecule,
                          POINT from,
                          double sph_step,
                          double max_radius,
                          dvector **fp);


/* Calculate the dihedral angle
 * for a given couple of atom (a1, a2) connected
 * by two common atoms (c1,c2):
 * a1-c1-c2-a2
 */
void Get3DDihedralAngle(MOLECULE molecule,
                       int atom_id1,
                       int atom_id2,
                       double *dh_angle,
                       uivector **dhids);

/* Calculate the dihedral angle
* for a given couple of atom (a1,a2) connected
* by a common atom (c)
* a1-c-a2
*/
void Get3DAngle(MOLECULE molecule,
               int atom_id1,
               int atom_id2,
               double *angle,
               uivector **aid);
#endif
