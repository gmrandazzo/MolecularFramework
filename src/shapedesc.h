/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef SHAPEDESC_H
#define SHAPEDESC_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <scientific.h>

#include "molecule.h"
#include "periodic_table.h"

/*
 * Shape point data structure
 */
typedef struct SHAPEPNT
{
  matrix *surf;
  matrix *vol;
} SHAPEPNT;

/* Calculate the shape points
 * - Volume points
 * - Surface points
 */
void getShapePoints(MOLECULE molecule, SHAPEPNT **pnt, enum RADIUS_TYPE rtype);

/*
 * Delete the shape points
 */
void deleteShapePoints(SHAPEPNT **pnt);

/*
 * Write the surface point in a mol2 file
 */
void WriteMol2SurfPoints(SHAPEPNT *pnt, char *outfile);

/*
 * Write the volume point in a mol2 file
 */
void WriteMol2VolPoints(SHAPEPNT *pnt, char *outfile);

/*
 * Calculate Van der Waals Molecular Volume
 */
void GetVDWMolVolSurf(MOLECULE molecule, double *volume, double *surface, double stepsize);

#endif
