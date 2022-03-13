/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include <stdio.h>

#include "molecule.h"
#include "periodic_table.h"

typedef struct {
  char *notes;
  char *atype_hash;
  char *atom_name;
  int hybridtype;
  double charge;
  double vdw_R_eq;
  double vdw_Epsilon;
  double hb_Umoment;
  double hb_Hcharge;
  size_t freq;
} Params;

typedef struct
{
  Params *params;
  size_t size;
} ForceField;

/*
 *  Load forcefield parameters from txt file
 */
void LoadFFParams(const char *ffparms, ForceField *ffield);

/*
 * Transfer  the mol2 molecule atom charges to the forcefield
 * and calculate the dipole in OH, NH, SH if present.
 */
void FFParamsAssigner(MOLECULE molecule, ForceField *ffield);

/*
 * Save the forcefield parameters to an output file
 */
void SaveFFParams(ForceField ffield, const char *outfile);

/*
 * Remove forcefield parameters;
 */
void DeleteForceField(ForceField *ffield);

/*
 * Search atom type in force field using atype_hash
 * N.B.: atype_hash is produced when import a molecule
 */
int atomtypeSearch(ForceField ffield, char *atype_hash, char *atom_name, int hybridtype);

/*
 * Assign force field parameters in molecule for future calculations
 */
void AssignParams(MOLECULE *molecule, ForceField ffield, enum RADIUS_TYPE rtype, int formal_charge);

void PrintFFParams(ForceField ff);

#endif
