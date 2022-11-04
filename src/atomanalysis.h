/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */
#ifndef ATOMANALYSIS_H
#define ATOMANALYSIS_H
#include <stdio.h>
#include <stdlib.h>
#include <scientific.h>

#include "molecule.h"

/*
 * Analyse molecule atom by atom assigning ATOMINFO parameters
 */
void AtomAnalyzer(MOLECULE *molecule, size_t fp_bond_max);

/*
 * Analyze atom connectivity and assign atom type
 */
void AtomTypeAnalyzer(MOLECULE *molecule);

/*
 * Generate the binary adjacence matrix
 * 1 means there is a connection
 * 0 means there is no connection
 */
void AdjMatGen(MOLECULE *molecule, matrix *adjmx);


/*
 * Generate the bond lenght adjacence matrix
 * If there is a bond between two atoms the value is > 0 and
 * represent the bond length.
 */
void BondLenghtAdjMatGen(MOLECULE *molecule, matrix *adjmx);

/*
 * Generate the coloured  adjacence matrix
 * Four kind of bond type are possible
 * 1, 2, 3, ar
 */
void BondColorAdjMatGen(MOLECULE *molecule, matrix *adjmx);

/*
*/
void Kekulize(MOLECULE *molecule);

#endif
