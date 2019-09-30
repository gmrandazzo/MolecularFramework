/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */
#include "molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Import some atom parameters from RDKit
 */
void rdkitMol2Read(const char *filename, MOLECULE *molecule);

/*
 * Search for rings
 */
void ringSearch(MOLECULE *molecule);

/*
 * Calculate the atom depth starting from a source atom.
 */
void getAtomDepth(MOLECULE *molecule, int src, int depth_max, int **atom_depth_);

/*
 * Convert atom to fingerprint for atomtype analysis
 * The fingerprint is function of atom types and bond length wich is controlled with depth_max.
 * Actually one bond length accounts for 15 atoms. See "GetFPPosition" function.
 */
void atomToFingerPrint(MOLECULE *molecule, int depth_max);

//void Align3DConformations(MOLECULE *from, MOLECULE *to);

#ifdef __cplusplus
} // extern "C"
#endif
