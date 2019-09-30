/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <math.h>
#include "misc.h"
#include "atomdesc.h"
#include "periodic_table.h"


void GetSandersonEneg(MOLECULE molecule, double *eneg)
{

}

void DipoleMoment(MOLECULE molecule, double *u, double *ux, double *uy, double *uz)
{
  int i;
  double mux, muy, muz;
  mux = muy = muz = 0.f;

  /*
   * the dipole moment is calculated trow these two formula:
   *
   * u(x) = coord_x * atom_charge
   *
   * mu = sqrt(u(x)^2 + u(y)^2 + u(z)^2)
   */
  for(i = 0; i < molecule.n_atoms; i++){
    mux += molecule.atoms[i].coord.x*molecule.atoms[i].charge;
    muy += molecule.atoms[i].coord.y*molecule.atoms[i].charge;
    muz += molecule.atoms[i].coord.z*molecule.atoms[i].charge;
  }
/*mux, muy and muz are expressed in e * Angstrom.
  * We convert these value in  Coulomb * metro (C * m)
  * e Angstrom * 1.6021765314E-19 (C) * 1.0E-10(A) = C*m
  * C*m is converted in debye dividing by 3.336·10^(-30) C·m
  * because 1D = 3.336·10^(-30) C·m
  */
  mux *= (1.6021765314e-29); mux /= 3.33564095198152e-30;
  muy *= (1.6021765314e-29); muy /= 3.33564095198152e-30;
  muz *= (1.6021765314e-29); muz /= 3.33564095198152e-30;

  /*From here dipole moments are in debye*/
  if(ux != NULL)
    (*ux) = mux;

  if(uy != NULL)
    (*uy) = muy;

  if(uz != NULL)
    (*uz) = muz;

  if(u != NULL)
    (*u) = sqrt(square(mux)+square(muy)+square(muz));
}

double GetTotalCharge(MOLECULE molecule)
{
  size_t i;
  double ch = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
    ch +=molecule.atoms[i].charge;
  }
  return ch;
}

/*
 *  If in a bond two atom have different electronegativity and one is bigger
 *  than the other, the most electronegative atom will get a
 *  delta - partial charge.
 * i.e.: C-H; eneg(C)=2.5 eneg(H)=2.1 2.5>2.1 then C(delta-), H(delta+)
 */
void GetElectronSpecie(MOLECULE molecule, int *pch_desc)
{
  size_t i, j;
  double eneg_centre = 0.f, eneg_neigh;
  for(i = 0; i < molecule.n_atoms; i++){
    eneg_centre = getEnegativity(molecule.atoms[i].asymbl);
    int eneg_neigh_big_than_centre = 0;
    for(j = 0; j < molecule.n_bonds; j++){
      if(molecule.bonds[j].origin_atom_id == i){
        eneg_neigh = getEnegativity(molecule.atoms[molecule.bonds[j].target_atom_id].asymbl);
        if(eneg_centre < eneg_neigh){
          eneg_neigh_big_than_centre = 1;
          break;
        }
        else{
          continue;
        }
      }
      else if(molecule.bonds[j].target_atom_id == i){
        eneg_neigh = getEnegativity(molecule.atoms[molecule.bonds[j].origin_atom_id].asymbl);
        if(eneg_centre < eneg_neigh){
          eneg_neigh_big_than_centre = 1;
          break;
        }
        else{
          continue;
        }
      }
      else{
        continue;
      }
    }
    if(eneg_neigh_big_than_centre == 1)
      pch_desc[i]= +1; /*delta+*/
    else
      pch_desc[i]= -1; /*delta-*/
  }
}


void GetInductiveEffects(MOLECULE molecule, int *pch_desc)
{
  size_t i, j;
  double eneg_centre = 0.f, eneg_neigh;
  for(i = 0; i < molecule.n_atoms; i++){
    eneg_centre = getEnegativity(molecule.atoms[i].asymbl);
    int eneg_neigh_big_than_centre = 0;
    for(j = 0; j < molecule.n_bonds; j++){
      if(molecule.bonds[j].origin_atom_id == i){
        eneg_neigh = getEnegativity(molecule.atoms[molecule.bonds[j].target_atom_id].asymbl);
        if(eneg_centre < eneg_neigh){
          eneg_neigh_big_than_centre = 1;
          break;
        }
        else{
          continue;
        }
      }
      else if(molecule.bonds[j].target_atom_id == i){
        eneg_neigh = getEnegativity(molecule.atoms[molecule.bonds[j].origin_atom_id].asymbl);
        if(eneg_centre < eneg_neigh){
          eneg_neigh_big_than_centre = 1;
          break;
        }
        else{
          continue;
        }
      }
      else{
        continue;
      }
    }
    if(eneg_neigh_big_than_centre == 1)
      pch_desc[i]= +1; /*delta+*/
    else
      pch_desc[i]= -1; /*delta-*/
  }
}
