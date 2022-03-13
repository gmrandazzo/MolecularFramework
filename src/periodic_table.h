/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

/*
 * Periodic table
 * Updated at 15 Jan 2014
 */

#include <stdio.h>
#include <string.h>
#include "misc.h"


/*
 * Enumerate the different radius type
 */
enum RADIUS_TYPE {
 vanderwaals,
 covalent
};

/*
 * Getting access to the periodic table information
 */
int RelAtomicMassTableSize();
void getRelAtomicMass(int i, char **atom, double *relams);

/*
 * Get covalent radius
 */
int CovRadiusTableSize();
void getCovRadius(int i, char **atom, double *cradius);
double getCovRadiusfromAtomName(char *_atom);

/* Van der Waals Radii
 * used to calculate the molecular volume
 */
double getVanDerWaalsRadiifromAtomName(char *_atom);

/*
 * Van der Waals Parameters from UFF Force Field
 * D_ii is expressed in expressed in kcal/mol
 * x_ii is expressed in Angstrom
 */
void getVanDerWaalsParams(char *_atom, double *D_ii, double *x_ii);

/*
 * Get First ionization potential for a given atom
 */
double getFirstIonizationPotential(char *_atom);

/*
 * Static scalar dipole polarizabilities (in atomic units) for neutral atoms
 * doi: 10.1080/00268976.2018.1535143
 */
double getStaticScalarDipolePolarizability(char *_atom);

/* get the atom number from atom name */
int getAtomNumberfromAtomName(char *_atom);


/* relative atomic mass
 * used to calculate the molecular weigth
 */
double getAtomWeightfromAtomName(char *_atom);


/* isotopic mass
 * used to calculate the exact mol weight
 */
double getExactAtomWeightfromAtomName(char *_atom);


/* Atom Lonepairs
 * This function will return to you the atom lone pair
 * 1: one lonepairs
 * 2: two lonepairs
 * 3: three lonepairs
 */
int getAtomLonepairs(char *_atom, int connectivity, int bonds);


/*
 * Get atom valence electrons
 */
int getAtomValenceElectrons(char *_atom);


/* Generic Atom Property dictionary read/search */
typedef struct{
  char type[16];
  double property;
} AProp;

typedef struct{
  size_t natoms;
  AProp *atoms;
} AtomsProperty;

void ReadAtomProperties(char *faprops, AtomsProperty **lst);
void DeleteAtomProperties(AtomsProperty **lst);
double getGenericAProperty(char *_type, AtomsProperty *lst);

/*
 * Get isotopes from atom name.
 */

/*
 * Isotope data structure definition
 */
#define MAX_ISOTEPES 9

typedef struct{
  char *name;
  int n_isotopes;
  double exactmass[MAX_ISOTEPES];
  double abundance[MAX_ISOTEPES];
} Isotope;

Isotope getIsotopefromAtomName(char *_atom);


double getEnegativity(const char *type);


/* return 0 if is not an atom and 1 if is an atom */
int isatom(char *_atom);


#endif
