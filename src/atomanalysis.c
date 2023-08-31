/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <math.h>
#include <scientific.h>
#include <unistd.h>

#include "atomanalysis.h"
#include "cppwrap.h"
#include "geomdesc.h"
#include "misc.h"
#include "periodic_table.h"

/*typedef struct ATYPES {
  MOLECULE *atypes
} ATYPES;

void ReadATypes(char *fatypes, ATYPES **atypes)
{
  (*atypes) = malloc(sizeof(ATYPES));
  (*atypes)->atypes = malloc(sizeof(MOLECULE)*35)
  NewEmptyMolecule((*atypes)->atypes[0], 2, 1) //sp2 C in C=O, C=S

  strcpy((*atypes)->atypes[0]->atoms[0].name, "C");
  strcpy((*atypes)->atypes[0]->atoms[1].name, "O");
  (*atypes)->atypes[0]->bonds[0].origin_atom_id = 0;
  (*atypes)->atypes[0]->bonds[0].target_atom_id = 1;
  strcpy((*atypes)->atypes[0]->bonds[0].type, "2");

}*/

/*
struct ATOMCONLONEPAIRS{
  int connectivity;
  int lonepairs;
  bool exception_rule;
  string atom_type;
};

const ATOMCONLONEPAIRS atomconlonepairs[] = {
  { 2, 2, true, "O"},
  { 1, 2, true, "O"},
  { 3, 1, true, "N"}, // Amines
  { 2, 1, true, "N"},// Exaple: R-N=O Group
  { 1, 1, true, "N"}, // Example: Nitrile group -CN
  { 2, 2, true, "S"},
  { 3, 1, false, "S"}, // no exception rule because the lonepairs is on d orbital, so no interaction could be possible to pi systems.
  { 3, 1, true, "P"}, // WARNING... maybe same to previous sulfur atom.
  { 1, 3, true, "F"},
  { 1, 3, true, "Cl"},
  { 1, 3, true, "Br"},
  { 1, 3, true, "I"}
};

int AtomHaveLonepairsException(MOLECULE *molecule, size_t aid){
  int i;
  const int sz =  sizeof(atomconlonepairs) / sizeof(atomconlonepairs[0]);

  for( i=0; i<sz; i++ ){
    if(atomconlonepairs[i].atom_type.compare(molecule->atoms[aid].type) == 0  &&
       atomconlonepairs[i].connectivity == molecule->atoms[aid].ainfo.connectivity){ // if the hatom have the lonepairs
      if(atomconlonepairs[i].exception_rule == true) // and if is subject to exception rule
        return 0;
      else
        return 1;
    }
  }
  return 1;
}*/

void AtomAnalyzer(MOLECULE *molecule, size_t fp_bond_max)
{
  /* Calcualte molecular connectivity */
  size_t i, j;
  for(i = 0; i < molecule->n_atoms; i++){
    /* Initializations of atominfo! */
    molecule->atoms[i].ainfo.n_Others = 0;
    molecule->atoms[i].ainfo.n_B = 0;
    molecule->atoms[i].ainfo.n_Cl = 0;
    molecule->atoms[i].ainfo.n_F = 0;
    molecule->atoms[i].ainfo.n_Br = 0;
    molecule->atoms[i].ainfo.n_I = 0;
    molecule->atoms[i].ainfo.n_P = 0;
    molecule->atoms[i].ainfo.n_S = 0;
    molecule->atoms[i].ainfo.n_N = 0;
    molecule->atoms[i].ainfo.n_O = 0;
    molecule->atoms[i].ainfo.n_C = 0;
    molecule->atoms[i].ainfo.n_H = 0;
    molecule->atoms[i].ainfo.aromatic = 0;

    /*Initializations of atom connectivity */
    molecule->atoms[i].ainfo.connectivity = 0;
    molecule->atoms[i].ainfo.bonds = 0;

    for(j = 0; j < molecule->n_bonds; j++){
      //printf("%lu %d %d %s\n", i, molecule->bonds[j].origin_atom_id , molecule->bonds[j].target_atom_id, molecule->bonds[j].type);
      if(molecule->bonds[j].origin_atom_id == i || molecule->bonds[j].target_atom_id == i){
        /* connectivity counter */
        molecule->atoms[i].ainfo.connectivity +=1;
        if(strcmp(molecule->bonds[j].type, "ar") != 0 ||
           strcmp(molecule->bonds[j].type, "am") != 0){
           molecule->atoms[i].ainfo.bonds += atoi(molecule->bonds[j].type);
        }

        /* hydrogen counter */
        size_t aid = 0;
        if(molecule->bonds[j].origin_atom_id == i){
          aid = molecule->bonds[j].target_atom_id;
        }
        else{ //else if(molecule->bonds[j].target_atom_id == i){
          aid = molecule->bonds[j].origin_atom_id;
        }

        char *atomsym = molecule->atoms[aid].asymbl;
        //strcpy(atomsym, molecule->atoms[aid].asymbl);

        if(strcmp(atomsym, "H") == 0){
          molecule->atoms[i].ainfo.n_H += 1;
        }
        else if(strcmp(atomsym, "C") == 0){
          molecule->atoms[i].ainfo.n_C += 1;
        }
        else if(strcmp(atomsym, "O") == 0){
          molecule->atoms[i].ainfo.n_O += 1;
        }
        else if(strcmp(atomsym, "N") == 0){
          molecule->atoms[i].ainfo.n_N += 1;
        }
        else if(strcmp(atomsym, "S") == 0){
          molecule->atoms[i].ainfo.n_S += 1;
        }
        else if(strcmp(atomsym, "P") == 0){
          molecule->atoms[i].ainfo.n_P += 1;
        }
        else if(strcmp(atomsym, "I") == 0){
          molecule->atoms[i].ainfo.n_I += 1;
        }
        else if(strcmp(atomsym, "Br") == 0){
          molecule->atoms[i].ainfo.n_Br += 1;
        }
        else if(strcmp(atomsym, "F") == 0){
          molecule->atoms[i].ainfo.n_F += 1;
        }
        else if(strcmp(atomsym, "Cl") == 0){
          molecule->atoms[i].ainfo.n_Cl += 1;
        }
        else if(strcmp(atomsym, "B") == 0){
          molecule->atoms[i].ainfo.n_B += 1;
        }
        else{
          /* WARNING! ATOM TYPE NOT SPECIFIED!*/
          molecule->atoms[i].ainfo.n_Others += 1;
        }
        //free(atomsym);
      }
      else{
        continue;
      }
    }

    /* Precalculation of lonepairs
     * At the end of the procedure you should recalculate the bonds for each atoms
     * and the lonepairs since aromaticity and kekule structure is generated.
     */
    molecule->atoms[i].ainfo.lonepairs = getAtomLonepairs(molecule->atoms[i].asymbl, molecule->atoms[i].ainfo.connectivity, molecule->atoms[i].ainfo.bonds);
    molecule->atoms[i].ainfo.valence_electrons = getAtomValenceElectrons(molecule->atoms[i].asymbl);

    /* Calculate atom hybridization
     * 0: unspecified; 1: s; 2: sp; 3: sp2; 4: sp3; 5:sp3d; 6: sp3d2; 7:unrecognized!
     */
     int t = molecule->atoms[i].ainfo.connectivity + molecule->atoms[i].ainfo.lonepairs;
     //printf("%s conn: %d lp: %d bonds: %d t: %d\n", molecule->atoms[i].name, molecule->atoms[i].ainfo.connectivity, molecule->atoms[i].ainfo.lonepairs, molecule->atoms[i].ainfo.bonds, t);
     switch(t){
      case 0:
        molecule->atoms[i].ainfo.hybrid = UNSPECIFIED;
        break;
      case 1:
        molecule->atoms[i].ainfo.hybrid = S;
        break;
      case 2:
        molecule->atoms[i].ainfo.hybrid = SP;
        break;
      case 3:
        molecule->atoms[i].ainfo.hybrid = SP2;
        break;
      case 4:
        molecule->atoms[i].ainfo.hybrid = SP3;
        break;
      case 5:
        molecule->atoms[i].ainfo.hybrid = SP3D;
        break;
      case 6:
        molecule->atoms[i].ainfo.hybrid = SP3D2;
        break;
      case OTHER:
        molecule->atoms[i].ainfo.hybrid = OTHER;
        break;
      default:
        molecule->atoms[i].ainfo.hybrid = UNSPECIFIED;
    }
  }

  /* Search for rings */
  ringSearch(molecule);

 /* Check for aromatic rings.
  * If ring with All SP2 atoms and follow the Huckel rule then each atom is AR
  */
  for(i = 0; i < molecule->n_rings; i++){
    /* calculate how many pi electrons we have from carbon atoms */
    for(j = 0; j < molecule->rings[i].size; j++){
      //printf("%s %d %d  ", molecule->atoms[molecule->rings[i].atoms[j]].asymbl, molecule->rings[i].atoms[j]+1, molecule->atoms[molecule->rings[i].atoms[j]].ainfo.hybrid);
      // Check if al carbons on ring are SP2
      if(strcmp(molecule->atoms[molecule->rings[i].atoms[j]].asymbl, "C") == 0)
        if(molecule->atoms[molecule->rings[i].atoms[j]].ainfo.hybrid == SP2 || molecule->atoms[molecule->rings[i].atoms[j]].ainfo.hybrid == AR)
          continue;
        else
          break;
      else
        continue;
    }


    if(j == molecule->rings[i].size){ // Maybe aromatic ring!

      size_t n_pi = 0;
      size_t n_lp = 0;

      for(j = 0; j < molecule->rings[i].size; j++){
        if(strcmp(molecule->atoms[molecule->rings[i].atoms[j]].asymbl, "C") == 0){ /* Pi carbon */
          if(molecule->atoms[molecule->rings[i].atoms[j]].ainfo.connectivity == 3)
            n_pi++;
        }
        else if(strcmp(molecule->atoms[molecule->rings[i].atoms[j]].asymbl, "N") == 0){
          //printf("N connectivity: %d\n", molecule->atoms[molecule->rings[i].atoms[j]].ainfo.connectivity );
          if(molecule->atoms[molecule->rings[i].atoms[j]].ainfo.connectivity == 2 ||
             molecule->atoms[molecule->rings[i].atoms[j]].ainfo.connectivity == 4){ /* Pi Nitrogen */
            n_pi++;
          }
          else if(molecule->atoms[molecule->rings[i].atoms[j]].ainfo.connectivity == 3){ /* lone painr nitrogen */
            n_lp++;
          }
          else{
            printf("WARNING! Nitrogen do not known!!\n");
          }
        }
        else if(strcmp(molecule->atoms[molecule->rings[i].atoms[j]].asymbl, "S") == 0){ /* one lonepair in ring */
          n_lp++;
        }
        else if(strcmp(molecule->atoms[molecule->rings[i].atoms[j]].asymbl, "O") == 0){ /* one lonepair in ring */
          n_lp++;
        }
        else{
          continue;
        }
      }


      if((n_pi%2) == 0){
        n_pi += 2*n_lp;
      }
      else{
        n_pi -=1;
        n_pi += 2*n_lp;
      }

      /* check if rings are planar
      matrix *coord;
      NewMatrix(&coord, molecule->rings[i].size, 3);
      for(j = 0; j < molecule->rings[i].size; j++){
        coord->data[j][0] = molecule->atoms[molecule->rings[i].atoms[j]].coord.x;
        coord->data[j][1] = molecule->atoms[molecule->rings[i].atoms[j]].coord.y;
        coord->data[j][2] = molecule->atoms[molecule->rings[i].atoms[j]].coord.z;
      }
      double planarity = CalcPlanarity(coord);
      printf("Planarity %f\n", planarity);
      DelMatrix(&coord); */

      //printf("Final n_pi: %zu\n", n_pi);
      /* Huckel rule
       * n_pi = 4N+2
       */
      double N = (n_pi-2.f)/4.f;

      if(j == molecule->rings[i].size &&  N == (int)N){
        //printf("Is aromatic\n");
        /* is aromatic! */
        for(j = 0; j < molecule->rings[i].size; j++){
          molecule->atoms[molecule->rings[i].atoms[j]].ainfo.hybrid = AR;
          molecule->atoms[molecule->rings[i].atoms[j]].ainfo.aromatic = 1;
        }

        /* Set all possible bonds belonging to the ring as ar */
        for(size_t k = 0; k < molecule->rings[i].size; k++){
          for(j = 0; j < molecule->n_bonds; j++){
            if((molecule->bonds[j].origin_atom_id == molecule->rings[i].atoms[k] &&
               molecule->bonds[j].target_atom_id == molecule->rings[i].atoms[k+1]) ||
               (molecule->bonds[j].target_atom_id == molecule->rings[i].atoms[k] &&
                molecule->bonds[j].origin_atom_id == molecule->rings[i].atoms[k+1])){
              strcpy(molecule->bonds[j].type, "ar");
              break;
            }
            else{
              continue;
            }
          }
        }
      }
      else{
        //printf("Not aromatic j: %zu rsize: %zu  nel: %zu N: %f\n", j,molecule->rings[i].size, n_pi, N);
        /*is not aromatic ring */
        continue;
      }
    }
  }

  /* set Hybridization for hydrogen */
  for(i = 0; i < molecule->n_bonds; i++){
    char *atomsym_orig = molecule->atoms[molecule->bonds[i].origin_atom_id].asymbl;
    char *atomsym_targ = molecule->atoms[molecule->bonds[i].target_atom_id].asymbl;
    if(strcmp(atomsym_orig, "H") == 0){
      molecule->atoms[molecule->bonds[i].origin_atom_id].ainfo.hybrid = 1;
      /*molecule->atoms[molecule->bonds[i].origin_atom_id].ainfo.hybrid = molecule->atoms[molecule->bonds[i].target_atom_id].ainfo.hybrid;*/
    }
    else if(strcmp(atomsym_targ, "H") == 0){
      molecule->atoms[molecule->bonds[i].target_atom_id].ainfo.hybrid = 1;
      /*molecule->atoms[molecule->bonds[i].target_atom_id].ainfo.hybrid = molecule->atoms[molecule->bonds[i].origin_atom_id].ainfo.hybrid;*/
    }
  }
  //PrintMolecule(*(molecule));

  //Kekulize(molecule);

  /*Load some parameters from RDKIT:
  - UFF or MMFF forcefield
  - Aromaticity
  - Stereoinformations
  */

  /* WARNING: FUNCTION NOT ACTIVE!
  rdkitMol2Read(mol2_file, molecule);
  */

  /* Calculate how many sp, sp2, sp3 or ar are bounded to each atom */
  for(i = 0; i < molecule->n_atoms; i++){
    molecule->atoms[i].ainfo.nsp = 0;
    molecule->atoms[i].ainfo.nsp2 = 0;
    molecule->atoms[i].ainfo.nsp3 = 0;
    molecule->atoms[i].ainfo.nar = 0;

    for(j = 0; j < molecule->n_bonds; j++){
      if(molecule->bonds[j].origin_atom_id == i || molecule->bonds[j].target_atom_id == i){
        size_t aid = 0;
        if(molecule->bonds[j].origin_atom_id == i){
          aid = molecule->bonds[j].target_atom_id;
        }
        else{ //else if(molecule->bonds[j].target_atom_id == i){
          aid = molecule->bonds[j].origin_atom_id;
        }

        if(molecule->atoms[aid].ainfo.hybrid == SP){
          molecule->atoms[i].ainfo.nsp += 1;
        }
        else if(molecule->atoms[aid].ainfo.hybrid == SP2){
          molecule->atoms[i].ainfo.nsp2 += 1;
        }
        else if(molecule->atoms[aid].ainfo.hybrid == SP3){
          molecule->atoms[i].ainfo.nsp3 += 1;
        }
        else if(molecule->atoms[aid].ainfo.hybrid == AR){
          molecule->atoms[i].ainfo.nar += 1;
        }
        else{ // else if(molecule->atoms[aid].hybrid == 4){
          continue;
        }
      }
    }
  }

  /*recalculate lonepairs and bonds for each atom since aromatic were detected!!*/

  /*
   * Convert atoms to a fingerprint to define Forcefield atom types
   * WARNING fp_bond_max
   * THIS IS AN IMPORTANT PARAMETER TO GENERATE THE ATOMTYPE HASH!
   * 1 means one bond distance, 2 two bond distance, etc...
   * N.B.: this value increase with the precision!
   */
  atomToFingerPrint(molecule, fp_bond_max);
}

void AdjMatGen(MOLECULE *molecule, matrix *adjmx)
{
  ResizeMatrix(adjmx, molecule->n_atoms, molecule->n_atoms);
  int i, j;
  for(i = 0; i < molecule->n_atoms; i++){
    adjmx->data[i][i] = 0;
    for(j = 0; j < molecule->n_bonds; j++){
     if(molecule->bonds[j].origin_atom_id == i){
       adjmx->data[i][molecule->bonds[j].target_atom_id] = 1;
     }
     else if(molecule->bonds[j].target_atom_id == i){
       adjmx->data[i][molecule->bonds[j].origin_atom_id] = 1;
     }
     else
      continue;
    }
  }
}


void BondLenghtAdjMatGen(MOLECULE *molecule, matrix *adjmx)
{
  ResizeMatrix(adjmx, molecule->n_atoms, molecule->n_atoms);
  int i, j;
  for(i = 0; i < molecule->n_atoms; i++){
    adjmx->data[i][i] = 0;
    for(j = 0; j < molecule->n_bonds; j++){
      if(molecule->bonds[j].origin_atom_id == i){
        int tid = molecule->bonds[j].target_atom_id;
        adjmx->data[i][molecule->bonds[j].target_atom_id] = GetDistance_(molecule->atoms[i].coord, molecule->atoms[tid].coord);
      }
      else if(molecule->bonds[j].target_atom_id == i){
        int oid = molecule->bonds[j].target_atom_id;
        adjmx->data[i][molecule->bonds[j].origin_atom_id] = GetDistance_(molecule->atoms[i].coord, molecule->atoms[oid].coord);
      }
      else{
        continue;
      }
    }
  }
}

void BondColorAdjMatGen(MOLECULE *molecule, matrix *adjmx)
{
  ResizeMatrix(adjmx, molecule->n_atoms, molecule->n_atoms);
  int i, j;
  for(i = 0; i < molecule->n_atoms; i++){
    adjmx->data[i][i] = 0;
    for(j = 0; j < molecule->n_bonds; j++){

      double bond_color = 0.f;
      if(strcmp(molecule->bonds[j].type, "3") == 0){
        bond_color = 3.f;
      }
      else if(strcmp(molecule->bonds[j].type, "2") == 0){
        bond_color = 2.f;
      }
      else if(strcmp(molecule->bonds[j].type, "1") == 0){
        bond_color = 1.f;
      }
      else if(strcmp(molecule->bonds[j].type, "ar") == 0){
        bond_color = 4.f;
      }
      else{
        bond_color = 9999;
      }

      if(molecule->bonds[j].origin_atom_id == i){
        adjmx->data[i][molecule->bonds[j].target_atom_id] = bond_color;
      }
      else if(molecule->bonds[j].target_atom_id == i){
        adjmx->data[i][molecule->bonds[j].origin_atom_id] = bond_color;
      }
      else{
        continue;
      }
    }
  }
}

int ring_contain_ring(int *ring1, size_t r1, int *ring2, size_t r2)
{
  size_t i, j;
  size_t found = 0;
  for(i = 0; i < r1; i++){
    for(j = 0; j < r2; j++){
      if(ring1[i] == ring2[j]){
        found++;
        break;
      }
      else{
        continue;
      }
    }
  }

  if(r1 < r2){
    if(found == r1){
      //r1 contained in r2
      return 0;
    }
    else{
      return -1;
    }
  }
  else{
    if(found == r2){
      //r2 contained in r1
      return 1; // get r1
    }
    else{
      return -1;
    }
  }
}

int isaromaticring(MOLECULE *molecule, size_t rid)
{
  size_t i;
  for(i = 0; i < molecule->rings[rid].size; i++){
    if(molecule->atoms[molecule->rings[rid].atoms[i]].ainfo.hybrid == AR){
      continue;
    }
    else
      break;
  }

  printf("is aromatic ring? %zu %zu\n", i, molecule->rings[rid].size);
  if(i == molecule->rings[rid].size)
    return 1;
  else
    return 0;
}

void univoque_ar_rings(MOLECULE *molecule, RINGS **urings, size_t *n_urings)
{
  size_t i, j, k;
  uivector *skip_id;
  initUIVector(&skip_id);
  RINGS *rings = molecule->rings;
  size_t n_rings = molecule->n_rings;
  for(i = 0; i < n_rings; i++){
    if(isaromaticring(molecule, i) == 1){
      for(j = i+1; j < n_rings; j++){
        int r = ring_contain_ring(rings[i].atoms, rings[i].size, rings[j].atoms, rings[j].size);
        if(r == 0){ //r1 in r2 and r2 > r1! get r2 and skip r1!
          if(isaromaticring(molecule, j) == 1){ // is aromatic ring
            UIVectorAppend(skip_id, i);
          }
          else{
            continue;
          }
        }
        else if(r == 1){ // r2 in r1 and r1 < r2! get r1 and skip r2!
          //UIVectorAppend(skip_id, j);
          continue;
        }
        else{
          // maybe unique ring!
          continue;
        }
      }
    }
    else{
      //UIVectorAppend(skip_id, i);
      continue;
    }
  }

  (*n_urings) = n_rings - skip_id->size;
  (*urings) = malloc(sizeof(RINGS)*(*n_urings));
  k = 0;
  for(i = 0; i < n_rings; i++){
    if(UIVectorHasValue(skip_id, i) == 0){
      continue;
    }
    else{
      (*urings)[k].atoms = malloc(sizeof(int)*rings[i].size);
      (*urings)[k].size = rings[i].size;
      for(j = 0; j < rings[i].size; j++){
        (*urings)[k].atoms[j] = rings[i].atoms[j];
      }
      k++;
    }
  }
}


int atom_can_have_db(int atom_bonds, int atom_lp, int atom_vel, int current_priority)
{
  printf("atom can have db priority: %d vel: %d lp: %d  bonds: %d = , atom_vel-(atom_bonds+atom_lp) %d or atom_vel-(atom_bonds) %d\n", current_priority, atom_vel, atom_lp, atom_bonds, atom_vel-(atom_bonds+atom_lp), atom_vel-(atom_bonds));
  if(current_priority == 0){
    if(atom_vel-(atom_bonds+atom_lp) >= 1){
      printf("return a 1\n");
      return 1; /* Yes this atom can have a double bond */
    }
    else{
      printf("return b 0\n");
      return 0; /* This atom can not have a double bond */
    }
  }
  else{
    if(atom_lp >=1 && current_priority == 1){
      if(atom_vel-(atom_bonds) >= 1){
        printf("return c 1\n");
        return 1; /*this atom have a lonepair for a double bond.... */
      }
      else{
        printf("return d 0\n");
        return 0; /* This atom can not have a double bond */
      }
    }
    else{
      if(atom_vel-(atom_bonds+atom_lp) >= 1){
        printf("return e 1\n");
        return 1; /* Yes this atom can have a double bond */
      }
      else{
        printf("return f 0\n");
        return 0; /* This atom can not have a double bond */
      }
    }
  }
}

int check_valence_problems_(MOLECULE *molecule, uivector *ring_id, uivector *tmp_bonds, uivector *priority, int current_priority)
{
  /*
   * if the ring atom can have a double bond and
   * there is no neighbors that can have a double bond
   * then we made a mistake.
   */
  size_t j, k;
  for(j = 0; j < ring_id->size; j++){
    /**/printf("Examination of id: %zu atom: %s nbonds: %zu vel: %d  lp: %d\n",
           ring_id->data[j]+1,
           molecule->atoms[ring_id->data[j]].asymbl,
           tmp_bonds->data[j],
           molecule->atoms[ring_id->data[j]].ainfo.valence_electrons,
           molecule->atoms[ring_id->data[j]].ainfo.lonepairs);

    size_t this_atom_vel = molecule->atoms[ring_id->data[j]].ainfo.valence_electrons;
    size_t this_atom_lp = 2*molecule->atoms[ring_id->data[j]].ainfo.lonepairs;
    if(atom_can_have_db(tmp_bonds->data[j],
                        this_atom_lp,
                        this_atom_vel, current_priority) == 1){
      int neigh_db = 0; // 1 = yes ; 0 = no.
      /* check if neighbors can have a double bond */
      for(k = 0; k < molecule->n_bonds; k++){
        if(molecule->bonds[k].origin_atom_id == ring_id->data[j]){
          size_t neigh_atom_vel = molecule->atoms[molecule->bonds[k].target_atom_id].ainfo.valence_electrons;
          size_t neigh_atom_lp = 2*molecule->atoms[molecule->bonds[k].target_atom_id].ainfo.lonepairs;
          size_t neigh_bonds;
          int indx = UIVectorIndexOf(ring_id, molecule->bonds[k].target_atom_id);
          if(indx > -1){
            neigh_bonds = tmp_bonds->data[indx];
          }
          else{
            neigh_bonds = molecule->atoms[molecule->bonds[k].target_atom_id].ainfo.bonds;
          }

          if(atom_can_have_db(neigh_bonds,
                              neigh_atom_lp,
                              neigh_atom_vel, current_priority) == 1){
            neigh_db = 1;
          }
          else{
            continue;
          }
        }
        else if(molecule->bonds[k].target_atom_id == ring_id->data[j]){
          size_t neigh_atom_vel = molecule->atoms[molecule->bonds[k].origin_atom_id].ainfo.valence_electrons;
          size_t neigh_atom_lp = 2*molecule->atoms[molecule->bonds[k].origin_atom_id].ainfo.lonepairs;
          size_t neigh_bonds;
          int indx = UIVectorIndexOf(ring_id, molecule->bonds[k].origin_atom_id);
          if(indx > -1){
            neigh_bonds = tmp_bonds->data[indx];
          }
          else{
            neigh_bonds = molecule->atoms[molecule->bonds[k].origin_atom_id].ainfo.bonds;
          }

          if(atom_can_have_db(neigh_bonds,
                              neigh_atom_lp,
                              neigh_atom_vel,
                              current_priority) == 1){
            neigh_db = 1;
          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }

      /*if neighbors can have a double bound neigh_db == 1 continue else break!*/
      if(neigh_db == 1){
        continue;
      }
      else{
        printf("Fail in while checking: %zu\n", ring_id->data[j]+1);
        break;
      }
    }
    else{
      continue;
    }
  }

  if(j == ring_id->size){
    return 1;
  }
  else{
    return 0;
  }
}


int check_valence_problems(MOLECULE *molecule, uivector *ring_id, uivector *tmp_bonds, uivector *priority, int current_priority)
{
  /*
   * if the ring atom can have a double bond and
   * there is no neighbors that can have a double bond
   * then we made a mistake.
   */
  size_t j;
  for(j = 0; j < ring_id->size; j++){
    /**/printf("Examination of id: %zu atom: %s nbonds: %zu vel: %d  lp: %d\n",
           ring_id->data[j]+1,
           molecule->atoms[ring_id->data[j]].asymbl,
           tmp_bonds->data[j],
           molecule->atoms[ring_id->data[j]].ainfo.valence_electrons,
           molecule->atoms[ring_id->data[j]].ainfo.lonepairs);

    size_t this_atom_vel = molecule->atoms[ring_id->data[j]].ainfo.valence_electrons;
    size_t this_atom_lp = 2*molecule->atoms[ring_id->data[j]].ainfo.lonepairs;
    if(atom_can_have_db(tmp_bonds->data[j],
                        this_atom_lp,
                        this_atom_vel, current_priority) == 1){
      return 0;
    }
    else{
      continue;
    }
  }
  return 1;
}


void ApplyKekuleDB(MOLECULE *molecule, uivector *ring_id, uivector *db_id)
{
  printf("APPLYING FINAL DOUBLE BONDS\n");
  size_t j, k;
  /*set single and double bonds*/
  for(k = 0; k < molecule->n_bonds; k++){
    if(UIVectorHasValue(ring_id, molecule->bonds[k].origin_atom_id) == 0 &&
       UIVectorHasValue(ring_id, molecule->bonds[k].target_atom_id) == 0 &&
       strcmp(molecule->bonds[k].type, "ar") == 0){
      strcpy(molecule->bonds[k].type, "1");
    }
    else{
      continue;
    }
  }

  for(j = 0; j < db_id->size; j+=2){
    printf("%zu->%zu\n", db_id->data[j]+1, db_id->data[j+1]+1);
  }

  for(j = 0; j < db_id->size; j+=2){
    for(k = 0; k < molecule->n_bonds; k++){
      if((db_id->data[j] == molecule->bonds[k].origin_atom_id &&
          db_id->data[j+1] == molecule->bonds[k].target_atom_id) ||
         (db_id->data[j] == molecule->bonds[k].target_atom_id &&
          db_id->data[j+1] == molecule->bonds[k].origin_atom_id)){
        strcpy(molecule->bonds[k].type, "2");
        printf("Set double bond to: %d %d %s\n", molecule->bonds[k].origin_atom_id+1, molecule->bonds[k].target_atom_id+1, molecule->bonds[k].type);
        break;
      }
      else{
        continue;
      }
    }
  }
}

void Kekulize(MOLECULE *molecule)
{
  size_t i, j, k, it;
  uivector *ring_id;
  uivector *priority, *tmp_priority;
  uivector *db_id, *final_db_id;
  uivector *bonds, *tmp_bonds;

  size_t n_urings;
  RINGS *urings;
  univoque_ar_rings(molecule, &urings, &n_urings);

  i = 0; j = 0; k = 0; it = 0;
  /* Filtered urings  */
  printf("Filtered ring..\n");
  for(i = 0; i < n_urings; i++){
    for(j = 0; j < urings[i].size; j++){
      printf("%d ", urings[i].atoms[j]+1);
    }
    printf("\n");
  }


  for(i = 0; i < n_urings; i++){
    initUIVector(&final_db_id);
    //Is an aromatic ring...
    NewUIVector(&ring_id, urings[i].size);
    NewUIVector(&priority, urings[i].size);
    NewUIVector(&tmp_priority, urings[i].size);
    NewUIVector(&bonds, urings[i].size);
    NewUIVector(&tmp_bonds, urings[i].size);


    for(j = 0; j < urings[i].size; j++){
      ring_id->data[j] =  urings[i].atoms[j];
    }

    /* calculate the number of bonds for each atom */
    for(j = 0; j < ring_id->size; j++){
      for(k = 0; k < molecule->n_bonds; k++){
        if(molecule->bonds[k].origin_atom_id == ring_id->data[j] ||
           molecule->bonds[k].target_atom_id == ring_id->data[j]){
          if(strcmp(molecule->bonds[k].type, "ar") == 0){
            bonds->data[j]++;
            tmp_bonds->data[j]++;
          }
          else{
            bonds->data[j] += atoi(molecule->bonds[k].type);
            tmp_bonds->data[j] += atoi(molecule->bonds[k].type);
          }
        }
        else{
          continue;
        }
      }
    }

    /* Calculate the priority */
    for(j = 0; j < ring_id->size; j++){
      if(strcmp(molecule->atoms[ring_id->data[j]].asymbl, "N") == 0){
        if(bonds->data[j] == 3){
          priority->data[j] += 1;
        }
        else if(bonds->data[j] <= 2){
          priority->data[j] = 0;
        }
        else{
          priority->data[j] = 2;
        }
      }
      else if(strcmp(molecule->atoms[ring_id->data[j]].asymbl, "C") == 0){
        if(bonds->data[j] <= 3){
          priority->data[j] = 0;
        }
        else{
          priority->data[j] = 2;
        }
      }
      else{
        priority->data[j] = 2;
      }
    }

    /*priority copy*/
    for(j = 0; j < priority->size; j++){
      tmp_priority->data[j] = priority->data[j];
    }

    /**/printf("RING under examination:\n");
    for(j = 0; j < ring_id->size; j++){
      printf("%zu %zu %d\n", ring_id->data[j]+1, bonds->data[j], molecule->atoms[ring_id->data[j]].ainfo.hybrid);
    }
    printf("#############\n");



    it = 0;
    while(it < ring_id->size){
      initUIVector(&db_id);
      int prev_db_id = -1;
      for(size_t p = 0; p < 2; p++){
        for(j = 0; j < ring_id->size; j++){
          size_t this_atom_vel = molecule->atoms[ring_id->data[j]].ainfo.valence_electrons;
          size_t this_atom_lp = 2*molecule->atoms[ring_id->data[j]].ainfo.lonepairs;

          size_t n_indx;
          if(j == ring_id->size-1){
            n_indx = 0;
          }
          else{
            n_indx = j+1;
          }

          size_t next_atom_vel = molecule->atoms[ring_id->data[n_indx]].ainfo.valence_electrons;
          size_t next_atom_lp = 2*molecule->atoms[ring_id->data[n_indx]].ainfo.lonepairs;

          printf("Can we have a double bond using priority %zu??\n", p);
          printf("id: %zu priority %zu atom: %s nbonds: %zu vel: %d  lp: %d\n",
                 ring_id->data[j]+1,
                 tmp_priority->data[j],
                 molecule->atoms[ring_id->data[j]].asymbl,
                 tmp_bonds->data[j],
                 molecule->atoms[ring_id->data[j]].ainfo.valence_electrons,
                 molecule->atoms[ring_id->data[j]].ainfo.lonepairs);
         printf("id: %zu priority %zu atom: %s nbonds: %zu vel: %d  lp: %d\n",
                ring_id->data[n_indx]+1,
                tmp_priority->data[n_indx],
                molecule->atoms[ring_id->data[n_indx]].asymbl,
                tmp_bonds->data[n_indx],
                molecule->atoms[ring_id->data[n_indx]].ainfo.valence_electrons,
                molecule->atoms[ring_id->data[n_indx]].ainfo.lonepairs);

          if(atom_can_have_db(tmp_bonds->data[j],
                              this_atom_lp,
                              this_atom_vel, p) == 1 &&
            atom_can_have_db(tmp_bonds->data[n_indx],
                             next_atom_lp,
                             next_atom_vel, p) == 1 &&
            tmp_priority->data[j] == p &&
            tmp_priority->data[j] == tmp_priority->data[n_indx] &&
            UIVectorHasValue(db_id, ring_id->data[j]) == 1 &&
            UIVectorHasValue(db_id, ring_id->data[n_indx]) == 1
            && prev_db_id != ring_id->data[j]){
            //this couple of atom could have a double bond
            UIVectorAppend(db_id, ring_id->data[j]);
            UIVectorAppend(db_id, ring_id->data[n_indx]);
            tmp_bonds->data[j]+=1;
            tmp_bonds->data[n_indx]+=1;
            prev_db_id = ring_id->data[n_indx];
            tmp_priority->data[j] = tmp_priority->data[n_indx] = 2;
            /**/printf("%zu %zu YES!\n", ring_id->data[j]+1, ring_id->data[n_indx]+1);
            printf("uuuuuuuuuuuuuuu\n");
          }
          else{
            //single bond
            /**/printf("%zu %zu NO!\n", ring_id->data[j]+1, ring_id->data[n_indx]+1);
            printf("uuuuuuuuuuuuuuu\n");
            continue;
          }
        }

        printf("END J LOOP\n");

        if(check_valence_problems(molecule, ring_id, tmp_bonds, tmp_priority, p) == 1){
          break;
        }
        else{
          /* set previous assigne double boinds to = 2 */
          for(j = 0; j < ring_id->size; j++)
            if(tmp_priority->data[j] == 0)
              tmp_priority->data[j]+=1;

          printf("PRIORITY++\n");
        }
      }

      printf("$$$$$$$$$$$$$$$\n");

      /**/printf("Possible Double bonds:\n");
      for(j = 0; j < db_id->size; j++){
        printf("%zu ", db_id->data[j]+1);
      }
      printf("\n");
      printf("-------------\n");

      if(check_valence_problems(molecule, ring_id, tmp_bonds, tmp_priority, 0) == 1){
        /**/printf("KEKULE OK!\n");
        printf("Possible double bonds\n");
        for(j = 0; j < db_id->size; j++){
          printf("%zu ", db_id->data[j]+1);
        }
        printf("\n");
        //sleep(3);

        if(final_db_id->size == 0){
          for(j = 0; j < db_id->size; j++)
            UIVectorAppend(final_db_id, db_id->data[j]);
        }
        else{
          if(db_id->size < final_db_id->size){
            DelUIVector(&final_db_id);
            NewUIVector(&final_db_id, db_id->size);
            for(j = 0; j < db_id->size; j++)
              final_db_id->data[j] = db_id->data[j];
          }
          else{
            continue;
          }
        }
        //DelUIVector(&db_id);
        //break;
      }

      //else{
        /* valence problem
         * shuffle!
        */
        //printf("NO KEKULE OR WRONG init!!\n");

        size_t last_id = ring_id->data[ring_id->size-1];
        size_t last_nbond = bonds->data[bonds->size-1];
        size_t last_priority = priority->data[priority->size-1];
        for(j = ring_id->size-1; j > 0; j-=1){
          ring_id->data[j] = ring_id->data[j-1];
          bonds->data[j] = bonds->data[j-1];
          tmp_bonds->data[j] = bonds->data[j];
          priority->data[j] = priority->data[j-1];
          tmp_priority->data[j] = priority->data[j];
        }

        ring_id->data[0] = last_id;
        bonds->data[0] = last_nbond;
        tmp_bonds->data[0] = bonds->data[0];
        priority->data[0] = last_priority;
        tmp_priority->data[0] = priority->data[0];

        /**/
        printf("SHUFFLE\n");
        for(j = 0; j < ring_id->size; j++){
          printf("%zu ", ring_id->data[j]+1);
        }
        printf("\n");
        //sleep(1);

      //}
      DelUIVector(&db_id);
      it++;
    }

    ApplyKekuleDB(molecule, ring_id, final_db_id);

    DelUIVector(&ring_id);
    DelUIVector(&bonds);
    DelUIVector(&tmp_bonds);
    DelUIVector(&priority);
    DelUIVector(&tmp_priority);
    DelUIVector(&final_db_id);
  }

  /*Check for last ar bonds... this should be converted as single bonds.. */
  for(i = 0; i < molecule->n_bonds; i++){
    if(strcmp(molecule->bonds[i].type, "ar") == 0){
      strcpy(molecule->bonds[i].type, "1");
    }
    else{
      continue;
    }
  }

  for(i = 0; i < n_urings; i++){
    free(urings[i].atoms);
  }
  free(urings);
}

void Kekulize_(MOLECULE *molecule){
  size_t i;
  uivector *a;
  uivector *d;
  initUIVector(&a);
  initUIVector(&d);
  for(i = 0; i < molecule->n_atoms; i++){
    if(molecule->atoms[i].ainfo.aromatic == 1){
      UIVectorAppend(a, i);
      UIVectorAppend(d, 0);
    }
    else{
      continue;
    }
  }


  for(i = 0; i < molecule->n_bonds; i++){
    int indx = UIVectorIndexOf(a, molecule->bonds[i].origin_atom_id);
    if(indx == -1)
      indx = UIVectorIndexOf(a, molecule->bonds[i].target_atom_id);

    if(indx > -1){
      if(strcmp(molecule->bonds[i].type, "ar") == 0){
        a->data[indx] += 1;
      }
      else{
        a->data[indx] += atoi(molecule->bonds[i].type);
      }
    }
    else{
      continue;
    }
  }

  DelUIVector(&d);
  DelUIVector(&a);
}
