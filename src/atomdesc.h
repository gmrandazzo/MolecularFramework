/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef ATOMDESC_H
#define ATOMDESC_H

#include "molecule.h"

/*
 * Calculate the Sanderson electronegativity using the Sanderson's principle.
 * Principle of Electronegativity Equalization:
 * When two or more atoms initially different in electronegativity combine
 * chemically, they adjust to have the same intermediate electronegativity
 * within the compound. This intermediate electronegativity is given by
 * the geometric mean of the individual electronegativities of
 * the component atoms.
 * The geometric mean of n numbers is obtained by multiplying all
 * of the atomic electronegativitities together
 * and taking the nth root of the product.
 */
void GetSandersonEneg(MOLECULE molecule, double *eneg);


/* Dipole Moment Calculation
 *
 */
void DipoleMoment(MOLECULE molecule, double *u, double *ux, double *uy, double *uz);

/*
 * Calculate the total charge in the molecule
 */
double GetTotalCharge(MOLECULE molecule);

void GetElectronSpecie(MOLECULE molecule, int *pch_desc);

void GetInductiveEffects(MOLECULE molecule, int *pch_desc);


#endif
