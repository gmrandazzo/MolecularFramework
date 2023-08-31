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

#include "miscdesc.h"

#include "atomanalysis.h"
#include "molecule.h"
#include "periodic_table.h"
#include "misc.h"

extern inline int CountAtoms(MOLECULE molecule, const char *asymbl)
{
  int i, natoms = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    if(strcmp(molecule.atoms[i].asymbl, asymbl) == 0){
      natoms += 1;
    }
  }
  return natoms;
}

extern inline int CountAtomTypes(MOLECULE molecule, const char *atype)
{
  int i, natoms = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    if(strcmp(molecule.atoms[i].type, atype) == 0){
      natoms += 1;
    }
  }
  return natoms;
}

int GetNCarbonAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "C");
}

int GetNSP3CarbonAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "C.3");
}

int GetNSP2CarbonAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "C.2");
}

int GetNSPCarbonAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "C.1");
}

int GetNArCarbonAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "C.ar");
}

int GetNCatCarbonAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "C.cat");
}

int GetNHydrogenAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "H");
}

int GetNNitrogenAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "N");
}

int GetNQuaternaryNitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.4");
}

int GetNSP3NitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.3");
}

int GetNPlanarNitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.pl3");
}

int GetNSP2NitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.2");
}

int GetNSPNitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.1");
}

int GetNArNitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.ar");
}

int GetNamNitrogenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "N.am");
}

int GetNOxygenAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "O");
}

int GetNSP3OxygenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "O.3");
}

int GetNSP2OxygenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "O.2");
}

int GetNCO2OxygenAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "O.co2");
}

int GetNSulfurAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "S");
}

int GetNSP3SulfurAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "S.3");
}

int GetSOSulfurAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "S.O");
}

int GetNSO2SulfurAtoms(MOLECULE molecule)
{
  return CountAtomTypes(molecule, "S.O2");
}

int GetNPhosphorousAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "P");
}

int GetNFluorineAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "F");
}

int GetNChlorineAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "Cl");
}

int GetNBromineAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "Br");
}

int GetNIodineAtoms(MOLECULE molecule)
{
  return CountAtoms(molecule, "I");
}

int GetNNonHydrogenAtoms(MOLECULE molecule)
{
  int i, natoms_noh = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    if(strcmp(molecule.atoms[i].asymbl, "H") != 0){
      natoms_noh += 1;
    }
  }
  return natoms_noh;
}

int GetNOHBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    int oi = molecule.bonds[i].origin_atom_id;
    int ti = molecule.bonds[i].target_atom_id;
    if((strcmp(molecule.atoms[oi].asymbl, "O") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "H") == 0) ||
       (strcmp(molecule.atoms[oi].asymbl, "H") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "O") == 0)){
      n += 1;
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNNHBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    int oi = molecule.bonds[i].origin_atom_id;
    int ti = molecule.bonds[i].target_atom_id;
    if((strcmp(molecule.atoms[oi].asymbl, "N") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "H") == 0) ||
       (strcmp(molecule.atoms[oi].asymbl, "H") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "N") == 0)){
      n += 1;
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNSHBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    int oi = molecule.bonds[i].origin_atom_id;
    int ti = molecule.bonds[i].target_atom_id;
    if((strcmp(molecule.atoms[oi].asymbl, "S") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "H") == 0) ||
       (strcmp(molecule.atoms[oi].asymbl, "H") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "S") == 0)){
      n += 1;
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNPHBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    int oi = molecule.bonds[i].origin_atom_id;
    int ti = molecule.bonds[i].target_atom_id;
    if((strcmp(molecule.atoms[oi].asymbl, "P") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "H") == 0) ||
       (strcmp(molecule.atoms[oi].asymbl, "H") == 0 &&
        strcmp(molecule.atoms[ti].asymbl, "P") == 0)){
      n += 1;
    }
    else{
      continue;
    }
  }
  return n;
}

int GetCCDoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if(strcmp(molecule.atoms[oi].asymbl, "C") == 0 &&
         strcmp(molecule.atoms[ti].asymbl, "C") == 0){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetCODoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[oi].asymbl, "C") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "O") == 0) ||
         (strcmp(molecule.atoms[oi].asymbl, "O") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "C") == 0) ){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNNNDoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if(strcmp(molecule.atoms[oi].asymbl, "N") == 0 &&
         strcmp(molecule.atoms[ti].asymbl, "N") == 0){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNNODoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[oi].asymbl, "N") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "O") == 0) ||
         (strcmp(molecule.atoms[oi].asymbl, "O") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "N") == 0) ){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNSODoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[oi].asymbl, "S") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "O") == 0) ||
         (strcmp(molecule.atoms[oi].asymbl, "O") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "S") == 0) ){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNPODoubleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "2") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[oi].asymbl, "P") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "O") == 0) ||
         (strcmp(molecule.atoms[oi].asymbl, "O") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "P") == 0) ){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetCCTripleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "3") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if(strcmp(molecule.atoms[oi].asymbl, "C") == 0 &&
         strcmp(molecule.atoms[ti].asymbl, "C") == 0){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetNCNTripleBonds(MOLECULE molecule)
{
  int i, n;
  n = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "3") == 0){
      int oi = molecule.bonds[i].origin_atom_id;
      int ti = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[oi].asymbl, "C") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "N") == 0) ||
         (strcmp(molecule.atoms[oi].asymbl, "N") == 0 &&
          strcmp(molecule.atoms[ti].asymbl, "C") == 0) ){
        n += 1;
      }
      else{
        continue;
      }
    }
    else{
      continue;
    }
  }
  return n;
}

int GetRotableBonds(MOLECULE molecule)
{
  int i, nrbonds=0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(strcmp(molecule.bonds[i].type, "1") == 0){
      size_t af = molecule.bonds[i].origin_atom_id;
      size_t at = molecule.bonds[i].target_atom_id;
      if((strcmp(molecule.atoms[af].asymbl, "H") == 0 || strcmp(molecule.atoms[at].asymbl, "H") == 0) ||  /* Exclude hydrogens */
         (strcmp(molecule.atoms[af].type, "C.2") == 0 && strcmp(molecule.atoms[at].type, "N.am") == 0) || /* Exclude amide bonds */
         (strcmp(molecule.atoms[af].type, "N.am") == 0 && strcmp(molecule.atoms[at].type, "C.2") == 0)){ /* Exclude amide bonds */
        continue;
      }
      else{
        nrbonds++;
      }
    }
    else{
      continue;
    }
  }
  return nrbonds;
}

int GetSumValenceElectrons(MOLECULE molecule)
{
  int i, nvel = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    nvel += getAtomValenceElectrons(molecule.atoms[i].asymbl);
  }
  return nvel;
}

double GetSumFirstIonizationPotetial(MOLECULE molecule)
{
  int i;
  double fip = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    fip += getFirstIonizationPotential(molecule.atoms[i].asymbl);
  }
  return fip;
}

int bondHBDGroup(MOLECULE molecule, int h_id, int het_id)
{
  /* - -O-H
   * - -N-H
   * - -S-H
   * - -H-F
   * - -H-Cl
   * - -P-H
  */
  if((strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "O") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "N") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "S") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "P") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "F") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "Cl") == 0) ||
     (strcmp(molecule.atoms[h_id].asymbl, "H") == 0 && strcmp(molecule.atoms[het_id].asymbl, "Si") == 0)){/* DOI: 10.1039/C3CC46048G */
    return 1;
  }
  else{
    return 0;
  }
}

int atomHBAGroup(const char *asymbl)
{
  if(strcmp(asymbl, "O") == 0 ||
     strcmp(asymbl, "N") == 0 ||
     strcmp(asymbl, "S") == 0 ||
     strcmp(asymbl, "P") == 0 ||
     strcmp(asymbl, "F") == 0 ||
     strcmp(asymbl, "Cl") == 0
  ){
    return 1;
  }
  else{
    return 0;
  }
}

int GetFormalCharge(MOLECULE molecule)
{
  /* forma charge = N valence electron - (N lonepairs + N bonds) */
  size_t i, j;
  double n_vel = 0;
  double n_nonbonds_e = 0.f;
  double n_bonds = 0.f;
  double charge = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    n_vel = (double)getAtomValenceElectrons(molecule.atoms[i].asymbl);
    //n_nonbonds_e = 2*getAtomLonepairs(molecule.atoms[i].asymbl, molecule.atoms[i].ainfo.connectivity);
    /* calculate the bonded electrons */
    n_bonds = 0;
    for(j = 0; j < molecule.n_bonds; j++){
      if(molecule.bonds[j].origin_atom_id == i || molecule.bonds[j].target_atom_id == i){
        //printf("type: %s atoi: %d\n", molecule.bonds[j].type, atoi(molecule.bonds[j].type));
        if(strcmp(molecule.bonds[j].type, "ar") == 0){
          if(strcmp(molecule.atoms[i].asymbl, "C")){
            n_bonds += 1.5;
          }
          else if(strcmp(molecule.atoms[i].asymbl, "O") == 0){
            n_bonds += 1.;
          }
          else{
            n_bonds += 1.5;
          }
        }
        else if(strcmp(molecule.bonds[j].type, "am") == 0){
          n_bonds += 1;
        }
        else{
          n_bonds += (double)atoi(molecule.bonds[j].type);
        }
      }
      else{
        continue;
      }
    }

    /*printf("n_bonds float: %f\n", floor(n_bonds_));*/

    /* 8 rule with exceptions to calculate the non bonded electrons */
    if(strcmp(molecule.atoms[i].asymbl, "H") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Li") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Na") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "K") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Rb") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Cs") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Fr") == 0 ||

       strcmp(molecule.atoms[i].asymbl, "Be") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Mg") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Ca") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Sr") == 0 ||
       strcmp(molecule.atoms[i].asymbl, "Ba") == 0){
        n_nonbonds_e = 0.f;
    }
    else if(strcmp(molecule.atoms[i].asymbl, "S") == 0){
      if(FLOAT_EQ(n_bonds, 6, 1e-2))
        n_nonbonds_e = 0.f;
      //else if(n_bonds == 4)
      else if(FLOAT_EQ(n_bonds, 4, 1e-2))
        n_nonbonds_e = 2.f;
      else{
        n_nonbonds_e = 8.f-(2.f*n_bonds);
      }
    }
    else if(strcmp(molecule.atoms[i].asymbl, "N") == 0){
      /* Case of C=N(=O)-OH in wich N is + and =O is -*/
      if(FLOAT_EQ(n_bonds, 5, 1e-2)){
      //if(n_bonds == 5){
        n_nonbonds_e = -1.f;
      }
      /* Case of quaternary amonium*/
      else if(FLOAT_EQ(n_bonds, 5, 1e-2)){
      //else if(n_bonds == 4){
        n_nonbonds_e = 0.f;
      }
      else{
        n_nonbonds_e = 8.f-(2.f*n_bonds);
      }
    }
    else if(strcmp(molecule.atoms[i].asymbl, "P") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "As") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "Sb") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "Bi") == 0){
      if(FLOAT_EQ(n_bonds, 5, 1e-2)){
      //if(n_bonds == 5){
        n_nonbonds_e = 0.f;
      }
      else{
        n_nonbonds_e = 8.f-(2.f*n_bonds);
      }
    }
    else if(strcmp(molecule.atoms[i].asymbl, "B") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "Al") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "Ga") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "In") == 0 ||
            strcmp(molecule.atoms[i].asymbl, "Tl") == 0){
        n_nonbonds_e = 0.f;
    }
    else{
      n_nonbonds_e = 8.f-(2.f*n_bonds);
    }

    /*printf("%zu %s %s %f\n", i+1, molecule.atoms[i].asymbl, molecule.atoms[i].type, n_bonds);

    printf("%s %f\n",  molecule.atoms[i].asymbl, n_bonds);
    printf("n_bonds: %f\n", n_bonds);*/

    charge += (n_vel-(n_nonbonds_e+n_bonds));

    printf("%s %s ar: %d hybrid: %d charge: %f bonds: %f non_bond_e: %f vel: %f\n", molecule.atoms[i].asymbl, molecule.atoms[i].type, molecule.atoms[i].ainfo.aromatic, molecule.atoms[i].ainfo.hybrid, n_vel-(n_nonbonds_e+n_bonds),  n_bonds, n_nonbonds_e, n_vel);
    puts("-------------------");
  }
  return (int)charge;
}


int GetNonBondingElectrons(char *asymbl, int n_bonds)
{
  int n_nonbonds_e = 0;
  /* 8 rule with exceptions to calculate the non bonded electrons */
  if(strcmp(asymbl, "H") == 0 ||
     strcmp(asymbl, "Li") == 0 ||
     strcmp(asymbl, "Na") == 0 ||
     strcmp(asymbl, "K") == 0 ||
     strcmp(asymbl, "Rb") == 0 ||
     strcmp(asymbl, "Cs") == 0 ||
     strcmp(asymbl, "Fr") == 0 ||

     strcmp(asymbl, "Be") == 0 ||
     strcmp(asymbl, "Mg") == 0 ||
     strcmp(asymbl, "Ca") == 0 ||
     strcmp(asymbl, "Sr") == 0 ||
     strcmp(asymbl, "Ba") == 0){
      n_nonbonds_e = 0;
  }
  else if(strcmp(asymbl, "S") == 0){
    if(n_bonds == 6)
      n_nonbonds_e = 0;
    else if(n_bonds == 4)
      n_nonbonds_e = 2;
    else{
      n_nonbonds_e = 8-(2*n_bonds);
    }
  }
  /*else if(strcmp(asymbl, "O") == 0){
    if(n_bonds == 3){
      n_nonbonds_e = 1;
    }
    else{
      n_nonbonds_e = 2;
    }
  }*/
  else if(strcmp(asymbl, "N") == 0){
    /* Case of C=N(=O)-OH in wich N is + and =O is -*/
    if(n_bonds == 5){
      n_nonbonds_e = -1;
    }
    /* Case of quaternary amonium*/
    else if(n_bonds == 4){
      n_nonbonds_e = 0;
    }
    else{
      n_nonbonds_e = 8-(2*n_bonds);
    }
  }
  else if(strcmp(asymbl, "P") == 0 ||
          strcmp(asymbl, "As") == 0 ||
          strcmp(asymbl, "Sb") == 0 ||
          strcmp(asymbl, "Bi") == 0){
    if(n_bonds == 5){
      n_nonbonds_e = 0;
    }
    else{
      n_nonbonds_e = 8-(2*n_bonds);
    }
  }
  else if(strcmp(asymbl, "B") == 0 ||
          strcmp(asymbl, "Al") == 0 ||
          strcmp(asymbl, "Ga") == 0 ||
          strcmp(asymbl, "In") == 0 ||
          strcmp(asymbl, "Tl") == 0){
      n_nonbonds_e = 0;
  }
  else{
    n_nonbonds_e = 8-(2*n_bonds);
  }
  return n_nonbonds_e;
}

int CalcAtomFormalCharge(MOLECULE molecule, int aid)
{
  size_t i;
  int n_vel = (double)getAtomValenceElectrons(molecule.atoms[aid].asymbl);
  int n_bonds = 0;
  int non_bond_e = 0;
  /*Calculate the number of bonds in the atom "aid"*/
  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].origin_atom_id == aid || molecule.bonds[i].target_atom_id == aid){
      if(strcmp(molecule.bonds[i].type, "ar") == 0){
        n_bonds += 1.f;
      }
      else{
        n_bonds += (double)atoi(molecule.bonds[i].type);
      }
    }
    else{
      continue;
    }
  }

  /* Rule exception if the atom "aid" is aromatic and the
   * bonds detected here are < of the valence electrons
   * then add +1 bond which represent the aromatic bond
   */
  if(strcmp(molecule.atoms[aid].asymbl, "C") == 0 &&
     molecule.atoms[aid].ainfo.aromatic == 1 &&
     n_bonds < n_vel){
    n_bonds += 1;
  }
  /*else if(strcmp(molecule.atoms[aid].asymbl, "O") == 0 &&
          molecule.atoms[aid].ainfo.aromatic == 1 &&
          n_bonds < n_vel){
    n_bonds -= 1;
  }*/

  /*The number of bonding electrons is 2*nbonds */
  non_bond_e = (double)GetNonBondingElectrons(molecule.atoms[aid].asymbl, n_bonds);

  printf("%d %s %s ar: %d hybrid: %d n_bonds: %d non_bond_e: %d vel: %d fc: %d\n", aid+1,
                                                                                   molecule.atoms[aid].asymbl,
                                                                                   molecule.atoms[aid].type,
                                                                                   molecule.atoms[aid].ainfo.aromatic,
                                                                                   molecule.atoms[aid].ainfo.hybrid,
                                                                                   n_bonds,
                                                                                   non_bond_e,
                                                                                   n_vel,
                                                                                   n_vel - non_bond_e - n_bonds);
  puts("-------------------");
  return n_vel - non_bond_e - n_bonds;
}

int GetNetCharge(MOLECULE molecule)
{
  size_t i;
  int charge = 0;
  for(i = 0; i < molecule.n_atoms; i++){
    charge += CalcAtomFormalCharge(molecule, i);
  }
  return charge;
}


double GetTotalPartialCharge(MOLECULE molecule)
{
  size_t i;
  double charge = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    charge += molecule.atoms[i].charge;
  }
  return charge;
}
