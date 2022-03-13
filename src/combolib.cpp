/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <iostream>
#include <vector>
#include "molecule.h"

using namespace std;

class ComboChemistry
{
public:
  ComboChemistry();
  ~ComboChemistry();
  std::vecror <MOLECULE> getComboMolecules();
  void setScaffold(std::string scfile);
  void setFramgent(std::string fragfile);
private:
  MOLECULE sc;
  std::vector<MOLECULE> frags, finalmol;
  std::vector <int> sc_du;

  void Combo(uint n, uint elements);
}

std::vecror <MOLECULE> ComboChemistry::getComboMolecules()
{
  std::vector combolst;
}

void ComboChemistry::setScaffold(std::string scfile)
{
  NewMolecule(&sc, scfile);
  sc_du.clear()
  for(int i = 0; i < sc->n_atoms; i++){
    if(strcmp(sc->atom[i]->type, "Du") == 0)
      sc_du.append(i);
    else
      continue;
  }
}

void ComboChemistry::setFramgent(std::string fragfile)
{
  frags.append(MOLECULE());
  NewMolecule(&frags.back(), fragfile);
}

ComboChemistry::ComboChemistry()
{

}

ComboChemistry::~ComboChemistry()
{
  for(int i = 0; i < frags.size(); i++)
    DelMolecule(&frags[i]);
  frags.clear();
  DelMolecule(&sc);
}

void AddMolecules(std::vector<MOLECULE> cmb)
{
  for(int i = 0; i < sc.n_bonds; i++){
    if(strcmp(sc.atoms[ sc.bonds[i].origin_atom_id-1], "Du")){

    }
    else if(strcmp(sc.atoms[sc.bonds[i].target_atom_id-1], "Du")){


    }
  }
}

void ComboChemistry::Combo()
{
  int n_dummy = sc_du.size(), n_fragments = frags.size();
  long long int i;
  int j;
  int d;
  long long int l = (long long int) powf((float) n_fragments, (float) n_dummy);
  std::vector<MOLECULE> comb;
  for (i = 0; i < n_dummy; i++)
    comb.append(0);

  for (i = 0; i < l; i++){
    d = 1;
    for (j = n_dummy-1; j > -1; j--){
      comb[j] = frags[(i/d) % n_fragments];
      d = d * n_fragments;
    }
    // comb contais the combination
    for(j = 0; j < n_dummy; j++){
      finalmol.append(MOLECULE());
      int n_atoms = ;
      int n_bonds = ;
      NewEmptyMolecule(&finalmol.back(), );
      finalmol.append();
      printf(" %c ", comb[j]);
    }
    printf("\n");
  }
}

//Combinazioni di lunghezza 2 di 4 elementi...
//n_cicli_for(2, 4);
