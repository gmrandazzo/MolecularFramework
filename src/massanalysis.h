/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef MASSANALYSIS_H
#define MASSANALYSIS_H

#include "molecule.h"

/*
 * Get the Nominal Mass
 * The nominal mass of a compound is its molecular weight
 * calculated using the atomic masses of constituent elements taken as integers.
 */

void GetNominalMW(MOLECULE molecule, int *mw);

/*
 * Get Molecular Weight
 */
void GetMW(MOLECULE molecule, double *mw);

/*
 * Get the Exact Molecular Weight
 */
void GetExactMW(MOLECULE molecule, double *exmw);

/*
 * Generate molecula formula
 */
void GetMolecularFormula(MOLECULE molecule, char **mformula);

/*
 * Generate molecular weight from molecula formula
 */
void GetMWfromMolecularFormula(char* mformula, double *mw);

/*
 * Calculate exact mass from molecular formula
 */
void GetExactMWfromMolecularFormula(char* mformula, double *exmw);

/*
 * Calculate the electronic chemical potential for solids state materials
 * using the Pauling Electronegativity
 * see doi:10.1016/j.chempr.2016.09.010
 */
double ElectronicChemicalPotentialFromFormula(char *formula);
double ElectronicChemicalPotentialFromMolecule(MOLECULE molecule);

/*
 * Calculate the isotopic distribution from molecula formula
 */
void GetIsotopicDistribution(char *mformula);

#endif
