/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef FIELDDESC_H
#define FIELDDESC_H
#include <stdio.h>
#include <stdlib.h>
#include <scientific.h>

#include "molecule.h"



void SphericalFieldAnalysis(MOLECULE molecule, matrix *field, double radius_step, matrix *spectra);

/* Single field descriptors */
double AvgFieldValue(matrix *field, double min, double max);
size_t FreqFieldValue(matrix *field, double min, double max);
double DipoleFieldValue(matrix *field, double min, double max);
double BaricenterFieldValue(MOLECULE molecule, matrix *field, double min, double max);

/* Double field descriptors */
double DiffFieldValue(matrix *field_a, matrix *field_b, double min, double max);
double RatioFieldValue(matrix *field_a, matrix *field_b, double min, double max);

#endif
