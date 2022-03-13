/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef MOLECULE_H
#define MOLECULE_H

#include <stdio.h>
#include <stdlib.h>

typedef struct{
  double x, y, z;
} POINT;

typedef enum hybridization {
  UNSPECIFIED = 0,
  S,
  SP,
  SP2,
  SP3,
  SP3D,
  SP3D2,
  AR,
  OTHER  //!< unrecognized hybridization
} hybridtype;

typedef struct{
  char *atype_hash;
  int valence_electrons;
  int lonepairs;
  int n_Others;
  int n_B;
  int n_Cl;
  int n_F;
  int n_Br;
  int n_I;
  int n_P;
  int n_S;
  int n_N;
  int n_O;
  int n_C;
  int n_H;
  int aromatic;
  hybridtype hybrid; /* see enum hybridization */
  int nsp;
  int nsp2;
  int nsp3;
  int nar;
  int stereo; /* 0: Nochiral 1: R, 2: S, 3: OTHER*/
  int cycle; /* 0: not cyclic; 1: cyclic */
  int connectivity;
  int bonds;
} ATOMINFO;

typedef struct {
  char name[256], type[256], asymbl[256], subst_name[256];
  int subst_id;
  POINT coord;
  ATOMINFO ainfo;
  double charge;
  double radius;
  int ffatompos; /*force field atom position*/
} ATOM;

typedef struct {
  char type[256];
  int origin_atom_id, target_atom_id;
} BOND;

typedef struct{
  int *atoms;
  size_t size;
} RINGS;

typedef struct {
  char molname[512];
  int n_atoms, n_bonds, n_rings;
  ATOM *atoms;
  BOND *bonds;
  RINGS *rings;
} MOLECULE;

void InitializeMolecule(MOLECULE *molecule);

void NewMOL2Molecule(MOLECULE *molecule , const char *filename);
void SaveMol2Molecule(MOLECULE molecule, const char *filename);
void NewSDFMolecule(MOLECULE *molecule, const char *filesdf);
void SaveSDFMolecule(MOLECULE *molecule, const char *filename);

void NewEmptyMolecule(MOLECULE *molecule, int n_atoms, int n_bonds);
void DelMolecule(MOLECULE *molecule);
void PrintMolecule(MOLECULE molecule);
void PrintMoleculeRings(MOLECULE molecule);
void PrintMoleculeAtomInfo(MOLECULE molecule);
void PrintMoleculeAtomFingerprint(MOLECULE molecule);

void SavePQRFile(MOLECULE molecule, const char *filename);

#endif
