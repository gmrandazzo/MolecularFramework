/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "kekulize.h"
#include "molecule.h"
#include <iostream>
#include <algorithm>
#include <vector>

int atom_can_have_db(int atom_bonds, int atom_lp, int atom_vel)
{
  //printf("atom can have db: %s vel: %d lp: %d  bonds: %d = %d\n", asymbl, atom_vel, atom_lp, atom_bonds, atom_vel-(atom_bonds+atom_lp));
  if(atom_vel-(atom_bonds+atom_lp) >= 1){
    return 1; /* Yes this atom can have a double bond */
  }
  else{
    return 0; /* This atom can not have a double bond */
  }
}

void Kekulize(MOLECULE *molecule)
{
  size_t i, j;
  std::vector<KATOM> atoms;

  for(i = 0; i < molecule->n_atoms; i++){
    if(molecule->atoms[i].ainfo.aromatic == 1){
      atoms.push_back(KATOM());
      atoms.back().aid = i;
      for(j = 0; j < molecule->n_bonds; j++){
        if(molecule->bonds[j].origin_atom_id == i){
          atoms.back().neighbours.push_back(molecule->bonds[j].target_atom_id);
        }
        else if(molecule->bonds[j].target_atom_id == i){
          atoms.back().neighbours.push_back(molecule->bonds[j].origin_atom_id);
        }
        else{
          continue;
        }
      }
    }
  }

  for(i = 0; i < atoms.size(); i++){
    int id = atoms[i].aid;

    if(atom_can_have_db(atoms[i].neighbours.size(),
                        molecule->atoms[id].ainfo.lonepairs,
                        molecule->atoms[id].ainfo.valence) == 1){
      for(j = 0; j < atoms[i].neighbours.size(); j++){
        size_t nid = atoms[i].neighbours[j];

        if(atom_can_have_db(atoms[i].neighbours.size(),
                            molecule->atoms[nid].ainfo.lonepairs,
                            molecule->atoms[nid].ainfo.valence) == 1){
      }
    }
    else{
      continue;
    }
  }
}
