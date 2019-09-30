/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef PROTANALYSIS_H
#define PROTANALYSIS_H
#include <scientific.h>

#include "molecule.h"

/* Find the structural water and return an atom id list*/
void findStructuralWater(MOLECULE molecule, uivector **atom_id);

/* Remove fromt the molecule the water except the structural ones
 * Algorithm:
 * Huggins et al
 * Protein Engineering, Design & Selection vol. 24 no. 10 pp. 777–789, 2011
 * DOI: 10.1093/protein/gzr036
 *
 */
void StripWater(MOLECULE inmolecule, MOLECULE *outmolecule);


#endif
