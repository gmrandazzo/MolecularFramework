/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef WFINGERPRINT_H
#define WFINGERPRINT_H

#include "molecule.h"

/*
 * Begin the Atom Fragment Model type.
 * This data structure is used to build fragment atom type and to calculate molecular propertyes from empirical assumption.
 */
typedef struct{
  char str[256];
  char atype[256];
  double aweight;
  int found;
} ATOMFRAG;

/*
 * Low level operation to get for each atom an atom fragment of type ATOMFRAG
 */
void GetAtomFrag(MOLECULE molecule, int aid, ATOMFRAG *afrag);

/*
 * This function is used to build an atom fragment table based on the actual imported molecule.
 */
void InitWeightTable(char *weightable, MOLECULE molecule);

/*
 * This function will calculate a molecular property by serching on the atom table the atom weight and summing these value
 *
 * Example usage:
 * argv[1] is a mol2 file
 * argv[2] is the table file.
 * MOLECULE molecule;
 * double property;
 * NewMolecule(&molecule, argv[1]);
 * CalcPropertFromSumAtomFrag(molecule, argv[2], &property);
 * printf("%s,%f\n", molecule.molname, property);
 * DelMolecule(&molecule);
 */
void CalcPropertFromSumAtomFrag(MOLECULE molecule, char *weightable, double *property);

/*
 * This function will return a fingerpint fp with size fpsize function of the weightable and each bit is sum of how many times
 * the atom fragment is found
 *
 * Example usage:
 * argv[1] is a mol2 file
 * argv[2] is the table file.
 *
 * double *fp = NULL;
 * int fpsize;
 * MOLECULE molecule;
 * NewMolecule(&molecule, argv[1]);
 * GetSumAtomFragFP(molecule, argv[2], &fp, &fpsize);
 *
 * for(int i = 0; i < fpsize; i++)
 *   printf("%f\n", fp[i])
 *
 * free(fp);
 *
 */
void GetSumAtomFragFP(MOLECULE molecule, char *weightable, double **fp, int *fpsize);

/*
 * This function will calculate a molecular property by serching on the atom table the atom weight and making the product of these value
 *
 * Example usage:
 * argv[1] is a mol2 file
 * argv[2] is the table file.
 * MOLECULE molecule;
 * double property;
 * NewMolecule(&molecule, argv[1]);
 * CalcPropertFromSumAtomFrag(molecule, argv[2], &property);
 * printf("%s,%f\n", molecule.molname, property);
 * DelMolecule(&molecule);
 */
void CalcPropertFromProductAtomFrag(MOLECULE molecule, char *strtable, double *property);

/*
 * This function will write for each molecule the corrispective atom fragment table
 */
void WriteAtomFragWeight(MOLECULE molecule, char *weightable, char *outcontribution);

#endif
