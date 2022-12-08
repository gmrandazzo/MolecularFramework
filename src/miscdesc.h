/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef MISCDESC_H
#define MISCDESC_H
#include <stdio.h>
#include <stdlib.h>

#include "molecule.h"

/*
 * Calculate the number of carbon atoms
 */
int GetNCarbonAtoms(MOLECULE molecule);

/*
 * N.B.: Depends on trypos atom type definition!
 */
int GetNSP3CarbonAtoms(MOLECULE molecule);

int GetNSP2CarbonAtoms(MOLECULE molecule);

int GetNSPCarbonAtoms(MOLECULE molecule);

int GetNArCarbonAtoms(MOLECULE molecule);

int GetNCatCarbonAtoms(MOLECULE molecule);

/*
 * Calculate the number of hydrogen atoms
 */
int GetNHydrogenAtoms(MOLECULE molecule);

/*  Calculate the number of nitrogen atoms
 *
 */
int GetNNitrogenAtoms(MOLECULE molecule);

/*
 * N.B.: Depends on trypos atom type definition!
 */
int GetNQuaternaryNitrogenAtoms(MOLECULE molecule);

int GetNSP3NitrogenAtoms(MOLECULE molecule);

int GetNPlanarNitrogenAtoms(MOLECULE molecule);

int GetNSP2NitrogenAtoms(MOLECULE molecule);

int GetNSPNitrogenAtoms(MOLECULE molecule);

int GetNArNitrogenAtoms(MOLECULE molecule);

int GetNamNitrogenAtoms(MOLECULE molecule);

/*
 * Calculate the number of oxygen atoms
 */
int GetNOxygenAtoms(MOLECULE molecule);

int GetNSP3OxygenAtoms(MOLECULE molecule);

int GetNSP2OxygenAtoms(MOLECULE molecule);

int GetNCO2OxygenAtoms(MOLECULE molecule);

/*
 * Calculate the number of sulfur atoms
 */
int GetNSulfurAtoms(MOLECULE molecule);

int GetNSP3SulfurAtoms(MOLECULE molecule);

int GetSOSulfurAtoms(MOLECULE molecule);

int GetNSO2SulfurAtoms(MOLECULE molecule);

/*
 * Calculate the number of phosphorous atoms
 */
int GetNPhosphorousAtoms(MOLECULE molecule);

/*
 * Calculate the number of alogen atoms
 */
 int GetNFluorineAtoms(MOLECULE molecule);

 int GetNChlorineAtoms(MOLECULE molecule);

 int GetNBromineAtoms(MOLECULE molecule);

 int GetNIodineAtoms(MOLECULE molecule);

/*
 * Calculate the number of non hydrogen atoms
 */
int GetNNonHydrogenAtoms(MOLECULE molecule);

/*
 * Calculate the number of -X-H bonds
 * with X = O, N, S, P
 */
int GetNOHBonds(MOLECULE molecule);

int GetNNHBonds(MOLECULE molecule);

int GetNSHBonds(MOLECULE molecule);

int GetNPHBonds(MOLECULE molecule);

/*
 * Calculate the number of X=X bonds
 * with X = C, N, S, P, O
 */

int GetCCDoubleBonds(MOLECULE molecule);

int GetCODoubleBonds(MOLECULE molecule);

int GetNNNDoubleBonds(MOLECULE molecule);

int GetNNODoubleBonds(MOLECULE molecule);

int GetNSODoubleBonds(MOLECULE molecule);

int GetNPODoubleBonds(MOLECULE molecule);

/*
 * Calculate the number of X#X bonds
 * with X = C, N
 */
int GetCCTripleBonds(MOLECULE molecule);

int GetNCNTripleBonds(MOLECULE molecule);


/*
 * Rotatable bondswere defined as:
 *  - any single bond, not ina ring,
 *  - bound to a nonterminal heavy (i.e., non-hydrogen) atom.
 * Excluded from the count were amide C-N bonds because oftheir high rotational energy barrier.
 *
 * Definition taken from
 * Veber et al J. Med. Chem.2002,45,2615-2623  doi: 10.1021/jm020017n
 *
 */
int GetRotableBonds(MOLECULE molecule);

/*
 * Calculate the sum of valence electron for the molecule
 */
int GetSumValenceElectrons(MOLECULE molecule);


/*
 * Calculate the sum of valence electron for the molecule
 */
double GetSumFirstIonizationPotetial(MOLECULE molecule);


/* Check if a bond is hydrogen bond donor group (HBD) */
int bondHBDGroup(MOLECULE molecule, int o_id, int t_id);

/* Check if an atom is hydrogen bond acceptor (HBA) */
int atomHBAGroup(const char *asymbl);

/* Calculate the molecular formal charge WRONG!!!*/
int GetFormalCharge(MOLECULE molecule);

int GetNetCharge(MOLECULE molecule);

double GetTotalPartialCharge(MOLECULE molecule);

#endif
