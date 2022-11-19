/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <scientific.h>

#include "atomanalysis.h"
#include "massanalysis.h"
#include "molecule.h"
#include "periodic_table.h"
#include "misc.h"

void GetNominalMW(MOLECULE molecule, int *mw)
{
  int i;
  (*mw) = 0;
  for(i = 0 ; i < molecule.n_atoms; i++){
    (*mw) += round(getAtomWeightfromAtomName(molecule.atoms[i].asymbl));
  }
}

void GetMW(MOLECULE molecule, double *mw)
{
  int i;
  (*mw) = 0.f;
  for(i = 0 ; i < molecule.n_atoms; i++){
    (*mw) += getAtomWeightfromAtomName(molecule.atoms[i].asymbl);
  }
}


void GetExactMW(MOLECULE molecule, double* exmw)
{
  int i;
  (*exmw) = 0.f;
  for(i = 0 ; i < molecule.n_atoms; i++){
    (*exmw) += getExactAtomWeightfromAtomName(molecule.atoms[i].asymbl);
  }
}

void GetMolecularFormula(MOLECULE molecule, char** mformula)
{

  int i, j, found;
  char *pch, int2char[256];
  (*mformula) = malloc(sizeof(char)*200);
  const int ntable = RelAtomicMassTableSize();
  char **aa = malloc(sizeof(char*)*ntable);
  for(i = 0; i < ntable; i++){
    aa[i] = malloc(sizeof(char)*3);
  }
  int *nn = malloc(sizeof(int)*ntable);

  strcpy(aa[0], "C");
  strcpy(aa[1], "H");
  strcpy(aa[2], "N");
  strcpy(aa[3], "P");
  strcpy(aa[4], "O");
  strcpy(aa[5], "S");
  strcpy(aa[6], "B");
  strcpy(aa[7], "F");
  strcpy(aa[8], "Cl");
  strcpy(aa[9], "Br");
  strcpy(aa[10], "I");
  strcpy(aa[11], "Si");

  j = 12;
  for(i = 0; i < ntable; i++){
    nn[i] = 0;
    char *atom_;
    getRelAtomicMass(i, &atom_, NULL);
    if(strcmp(atom_, aa[0]) != 0 &&
      strcmp(atom_, aa[1]) != 0 &&
      strcmp(atom_, aa[2]) != 0 &&
      strcmp(atom_, aa[3]) != 0 &&
      strcmp(atom_, aa[4]) != 0 &&
      strcmp(atom_, aa[5]) != 0 &&
      strcmp(atom_, aa[6]) != 0 &&
      strcmp(atom_, aa[7]) != 0 &&
      strcmp(atom_, aa[8]) != 0 &&
      strcmp(atom_, aa[9]) != 0 &&
      strcmp(atom_, aa[10]) != 0 &&
      strcmp(atom_, aa[11]) != 0){
      strcpy(aa[j], atom_);
      j++;
    }
    else{
      continue;
    }
  }

  for(i = 0; i < molecule.n_atoms; i++){
    pch = molecule.atoms[i].asymbl;
    found = -1;

    if(strcmp(pch, aa[0]) == 0){
      found = 0;
    }
    else if(strcmp(pch, aa[1]) == 0){
      found = 1;
    }
    else if(strcmp(pch, aa[2]) == 0){
      found = 2;
    }
    else if(strcmp(pch, aa[3]) == 0){
      found = 3;
    }
    else if(strcmp(pch, aa[4]) == 0){
      found = 4;
    }
    else if(strcmp(pch, aa[5]) == 0){
      found = 5;
    }
    else if(strcmp(pch, aa[6]) == 0){
      found = 6;
    }
    else if(strcmp(pch, aa[7]) == 0){
      found = 7;
    }
    else if(strcmp(pch, aa[8]) == 0){
      found = 8;
    }
    else if(strcmp(pch, aa[9]) == 0){
      found = 9;
    }
    else if(strcmp(pch, aa[10]) == 0){
      found = 10;
    }
    else if(strcmp(pch, aa[11]) == 0){
      found = 11;
    }
    else{
      for(j = 12; j < ntable; j++){
        if(strcmp(pch, aa[j]) == 0){
          found = j;
          break;
        }
        else{
          continue;
        }
      }
    }

    if(found == -1){
      fprintf(stderr,"Error! Uknown atom type: %s\n", molecule.atoms[i].type);
      return;
    }
    else{
      nn[found] +=1;
    }
  }

  if(nn[0] > 0){
    strcpy((*mformula), aa[0]);
    //sprintf(int2char, "%d", nn[0]); /*conversion integer to char */
    snprintf(int2char, sizeof(int2char), "%d", nn[0]);
    strcat((*mformula), int2char); /* appending the char */
  }
  else{
    strcpy((*mformula), "");
  }

  for(i = 1; i < ntable; i++){
    if(nn[i] > 0){
      strcat((*mformula), aa[i]);
      if(nn[i] > 1){
        //sprintf(int2char, "%d", nn[i]); /*conversion integer to char */
        snprintf(int2char, sizeof(int2char), "%d", nn[i]);
        strcat((*mformula), int2char); /* appending the char */
      }
    }
    else{
      continue;
    }
  }


  for(i = 0; i < ntable; i++){
    free(aa[i]);
  }
  free(aa);
  free(nn);
}

/*
 * Split the molecular formula into a list of string vector and int vector to describe
 * how many times the atom is counted.
 * N.B.: atom_names and n_atoms must be initialised.
 */
void ParseMolecularFormula(char* mformula, strvector **atom_names, ivector **n_atoms)
{
  int i, ait = 0, nit = 0;
  double na;
  char atom[100], n_atom[100];
  int previous_is_number = 1;

  const int sz = strlen(mformula);

  memset(atom, '\0', 100);
  memset(n_atom, '\0', 100);

  if((*atom_names) == NULL)
    initStrVector(atom_names);

  if((*n_atoms) == NULL)
    initIVector(n_atoms);

  for(i = 0; i < sz; i++){
    if(isalpha(mformula[i]) != 0){
      if(previous_is_number == 0){ /* Write and free the previous atom  name and add the new atom name*/
        /*Do calculation!*/
        atom[ait] = '\0';
        n_atom[nit] = '\0';
        na = atof(n_atom);
        StrVectorAppend((*atom_names), atom);
        if(FLOAT_EQ(na, 0, 1e-4))
          IVectorAppend((*n_atoms), 1);
        else
          IVectorAppend((*n_atoms), na);

        memset(atom, '\0', 100);
        memset(n_atom, '\0', 100);
        ait = 0;
        nit = 0;
        atom[ait] = mformula[i]; ait++;
        previous_is_number = 1;
      }
      else{ // append the atom name
        if(islower(mformula[i]) != 0){
          atom[ait] = mformula[i]; ait++;
        }
        else{
          if(isatom(atom) == 1){
            /*Do calculation!*/
            na = atof(n_atom);
            StrVectorAppend((*atom_names), atom);
            if(FLOAT_EQ(na, 0, 1e-4))
              IVectorAppend((*n_atoms), 1);
            else
              IVectorAppend((*n_atoms), na);


            memset(atom, '\0', 100);
            memset(n_atom, '\0', 100);
            ait = 0;
            nit = 0;
            atom[ait] = mformula[i]; ait++;
            previous_is_number = 1;
          }
          else{
            atom[ait] = mformula[i]; ait++;
          }
        }
      }
    }
    else if(isdigit(mformula[i]) != 0){
      previous_is_number = 0;
      n_atom[nit] = mformula[i]; nit++;
    }
    else{
      continue;
    }
  }

  na = atof(n_atom);
  StrVectorAppend((*atom_names), atom);
  if(FLOAT_EQ(na, 0, 1e-4))
    IVectorAppend((*n_atoms), 1);
  else
    IVectorAppend((*n_atoms), na);
}

void GetMWfromMolecularFormula(char* mformula, double *mw)
{
  int i = 0;
  strvector *atom_names;
  ivector *n_atoms;
  initStrVector(&atom_names);
  initIVector(&n_atoms);

  ParseMolecularFormula(mformula, &atom_names, &n_atoms);
  (*mw) = 0.f;
  for(i = 0; i < atom_names->size; i++){
    (*mw) += getAtomWeightfromAtomName(atom_names->data[i])*n_atoms->data[i];
  }
  DelStrVector(&atom_names);
  DelIVector(&n_atoms);
}

/*
void GetMWfromMolecularFormula(char* mformula, double *mw)
{
  int i, ait = 0, nit = 0;
  double na;
  char atom[100], n_atom[100];
  int previous_is_number = 1;

  const int sz = strlen(mformula);

  memset(atom, '\0', 100);
  memset(n_atom, '\0', 100);
  (*mw) = 0.f;
  for(i = 0; i < sz; i++){
    if(isalpha(mformula[i]) != 0){
      if(previous_is_number == 0){ // Write and free the previous atom  name and add the new atom name
        //Do calculation!
        atom[ait] = '\0';
        n_atom[nit] = '\0';
        na = atof(n_atom);
        if(FLOAT_EQ(na, 0, 1e-4))
          (*mw) += getAtomWeightfromAtomName(atom);
        else
          (*mw) += getAtomWeightfromAtomName(atom)*na;

        memset(atom, '\0', 100);
        memset(n_atom, '\0', 100);
        ait = 0;
        nit = 0;
        atom[ait] = mformula[i]; ait++;
        previous_is_number = 1;
      }
      else{ // append the atom name
        if(islower(mformula[i]) != 0){
          atom[ait] = mformula[i]; ait++;
        }
        else{
          if(isatom(atom) == 1){
            // Do calculation!
            na = atof(n_atom);
            if(FLOAT_EQ(na, 0, 1e-4))
              (*mw) += getAtomWeightfromAtomName(atom);
            else
              (*mw) += getAtomWeightfromAtomName(atom)*na;

            memset(atom, '\0', 100);
            memset(n_atom, '\0', 100);
            ait = 0;
            nit = 0;
            atom[ait] = mformula[i]; ait++;
            previous_is_number = 1;
          }
          else{
            atom[ait] = mformula[i]; ait++;
          }
        }
      }
    }
    else if(isdigit(mformula[i]) != 0){
      previous_is_number = 0;
      n_atom[nit] = mformula[i]; nit++;
    }
    else{P
      continue;
    }
  }

  na = atof(n_atom);
  if(FLOAT_EQ(na, 0, 1e-4))
    (*mw) += getAtomWeightfromAtomName(atom);
  else
    (*mw) += getAtomWeightfromAtomName(atom)*na;
}
*/

void GetExactMWfromMolecularFormula(char* mformula, double *mw)
{
  int i = 0;
  strvector *atom_names;
  ivector *n_atoms;
  initStrVector(&atom_names);
  initIVector(&n_atoms);

  ParseMolecularFormula(mformula, &atom_names, &n_atoms);
  (*mw) = 0.f;
  for(i = 0; i < atom_names->size; i++){
    (*mw) += getExactAtomWeightfromAtomName(atom_names->data[i])*n_atoms->data[i];
  }
  DelStrVector(&atom_names);
  DelIVector(&n_atoms);
}


double ElectronicChemicalPotentialFromFormula(char *formula)
{
  int i = 0;
  double X = 0.f, root = 0.f;
  strvector *atom_names;
  ivector *n_atoms;
  initStrVector(&atom_names);
  initIVector(&n_atoms);

  ParseMolecularFormula(formula, &atom_names, &n_atoms);
  for(i = 0; i < atom_names->size; i++){
    X += pow(getEnegativity(atom_names->data[i]),n_atoms->data[i]);
    root += (double)n_atoms->data[i];
  }
  DelStrVector(&atom_names);
  DelIVector(&n_atoms);
  return pow(X, 1.f/root);
}


double ElectronicChemicalPotentialFromMolecule(MOLECULE molecule)
{
  char *mform;
  GetMolecularFormula(molecule, &mform);
  double ret =  ElectronicChemicalPotentialFromFormula(mform);
  free(mform);
  return ret;
}

/*void GetExactMWfromMolecularFormula(char* mformula, double* exmw)
{
  int i, ait = 0, nit = 0;
  double na;
  char atom[100], n_atom[100];
  int previous_is_number = 1;

  const int sz = strlen(mformula);

  memset(atom, '\0', 100);
  memset(n_atom, '\0', 100);

  (*exmw) = 0.f;
  for(i = 0; i < sz; i++){
    if(isalpha(mformula[i]) != 0){
      if(previous_is_number == 0){ // Write and free the previous atom  name and add the new atom name

        // Do calculation!
        atom[ait] = '\0';
        n_atom[nit] = '\0';
        na = atof(n_atom);
        if(FLOAT_EQ(na, 0, 1e-4))
          (*exmw) += getExactAtomWeightfromAtomName(atom);
        else
          (*exmw) += getExactAtomWeightfromAtomName(atom)*na;

        memset(atom, '\0', 100);
        memset(n_atom, '\0', 100);
        ait = 0;
        nit = 0;
        atom[ait] = mformula[i]; ait++;
        previous_is_number = 1;
      }
      else{ // append the atom name
        if(islower(mformula[i]) != 0){
          atom[ait] = mformula[i]; ait++;
        }
        else{
          if(isatom(atom) == 1){
            // Do calculation!
            na = atof(n_atom);
            if(FLOAT_EQ(na, 0, 1e-4))
              (*exmw) += getExactAtomWeightfromAtomName(atom);
            else
              (*exmw) += getExactAtomWeightfromAtomName(atom)*na;

            memset(atom, '\0', 100);
            memset(n_atom, '\0', 100);
            ait = 0;
            nit = 0;
            atom[ait] = mformula[i]; ait++;
            previous_is_number = 1;
          }
          else{
            atom[ait] = mformula[i]; ait++;
          }
        }
      }
    }
    else if(isdigit(mformula[i]) != 0){
      previous_is_number = 0;
      n_atom[nit] = mformula[i]; nit++;
    }
    else{
      continue;
    }
  }

  na = atof(n_atom);
  if(FLOAT_EQ(na, 0, 1e-4))
    (*exmw) += getExactAtomWeightfromAtomName(atom);
  else
    (*exmw) += getExactAtomWeightfromAtomName(atom)*na;
}*/

void GetIsotopicDistribution(char* mformula)
{
  /*size_t i, j, k ,l, m;
  strvector *atom_names;
  ivector *n_atoms;
  initStrVector(&atom_names);
  initIVector(&n_atoms);

  ParseMolecularFormula(mformula, &atom_names, &n_atoms);

  double a_prev = 0.f, a;
  int k_prev, l_prev, m_prev;
  for(i = 0; i < atom_names->size; i++){
    Isotope iso = getIsotopefromAtomName(atom_names->data[i]);
    for(j = 0; j < iso.n_isotopes; j++){
      k = n_atoms->data[i];
      while(k != 0){
        l = n_atoms->data[i] - k;
        while(l != 0){
          m = n_atoms->data[i] - k - l;
          while(m != 0){
            if(a_prev > 0){
              a = a_prev * (Factorial(k_prev)*Factorial(l_prev)*Factorial(m_prev)) /(Factorial(k)*Factorial(l)*Factorial(m));
              //a *= iso;
            }
            else{
            }
          }
        }
        k--;
      }
    }
    a = a_prev*
    //(*mw) += getExactAtomWeightfromAtomName(atom_names->data[i])*n_atoms->data[i];
  }
  DelStrVector(&atom_names);
  DelIVector(&n_atoms);
  */
  return;
}
