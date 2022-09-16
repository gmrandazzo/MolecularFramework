/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "forcefield.h"
#include "molecule.h"
#include "atomanalysis.h"
#include "periodic_table.h"

#include "misc.h"
#include "miscdesc.h"
#include "atomdesc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>


void LoadFFParams(const char *ffparms, ForceField *ffield)
{
  int i, j;
  char l[2048], *tl;
  FILE *fffile;

  /* Allocate file in memory */
  if((fffile = fopen(ffparms, "r")) != NULL){
    i = 0;
    fseek(fffile, 0, SEEK_SET);
    while(fgets(l, LINE_MAX, fffile) != NULL){
      if(strstr(l, "#") != NULL){
        continue;
      }
      else{
        i++;
      }
    }

    ffield->params = malloc(sizeof(Params)*i);
    ffield->size = i;
    fclose(fffile);

    /* Charge the information in memory */
    fffile = fopen(ffparms, "r");

    i = 0;
    fseek(fffile, 0, SEEK_SET);
    while(fgets(l, LINE_MAX, fffile) != NULL){
      if (strstr(l, "#") != NULL){
        continue;
      }
      else{
        // Tokenize
        for(tl = strtok(Trim(l), " \t"), j = 0; tl; tl = strtok (NULL, " \t"), j++){
          if(j == 0){
            ffield->params[i].atype_hash = strdup(tl);
          }
          else if(j == 1){
            ffield->params[i].atom_name = strdup(tl);
          }
          else if(j == 2){
            ffield->params[i].hybridtype = atoi(tl);
          }
          else if(j == 3){
            ffield->params[i].charge = atof(tl);
          }
          else if(j == 4){
            ffield->params[i].vdw_R_eq = atof(tl);
          }
          else if(j == 5){
            ffield->params[i].vdw_Epsilon = atof(tl);
          }
          else if(j == 6){
            ffield->params[i].hb_Umoment = atof(tl);
          }
          else if(j == 7){
            ffield->params[i].hb_Hcharge = atof(tl);
          }
          else if(j == 8){
            ffield->params[i].freq = atoi(tl);
          }
          else if(j == 9){
            ffield->params[i].notes = strdup(tl);
          }
          else{
            continue;
          }
        }
        i++;
      }
    }
    fclose(fffile);
  }
  else{
    /* The file not exist but you can create a simple forcefile with 0 allocation */
    ffield->params = malloc(sizeof(Params)*0);
    ffield->size = 0;
  }
}

void FFHBParamsAssigner(MOLECULE molecule, int h_id, double *u, double *h_charge)
{
  MOLECULE umol;
  NewEmptyMolecule(&umol, 2, 0);
  size_t j;
  for(j = 0; j < molecule.n_bonds; j++){
    int oid = molecule.bonds[j].origin_atom_id;
    int tid = molecule.bonds[j].target_atom_id;
    if(oid == h_id){
      if(bondHBDGroup(molecule, oid, tid) == 1 || bondHBDGroup(molecule, tid, oid) == 1){
        umol.atoms[0].charge = molecule.atoms[oid].charge;
        umol.atoms[0].coord.x = molecule.atoms[oid].coord.x;
        umol.atoms[0].coord.y = molecule.atoms[oid].coord.y;
        umol.atoms[0].coord.z = molecule.atoms[oid].coord.z;

        umol.atoms[1].charge = molecule.atoms[tid].charge;
        umol.atoms[1].coord.x = molecule.atoms[tid].coord.x;
        umol.atoms[1].coord.y = molecule.atoms[tid].coord.y;
        umol.atoms[1].coord.z = molecule.atoms[tid].coord.z;
        DipoleMoment(umol, u, NULL, NULL, NULL);
        (*h_charge) = molecule.atoms[h_id].charge;
        break;
      }
    }
    else if(tid == h_id){
      if(bondHBDGroup(molecule, tid, oid) == 1 || bondHBDGroup(molecule, oid, tid) == 1){
        umol.atoms[0].charge = molecule.atoms[oid].charge;
        umol.atoms[0].coord.x = molecule.atoms[oid].coord.x;
        umol.atoms[0].coord.y = molecule.atoms[oid].coord.y;
        umol.atoms[0].coord.z = molecule.atoms[oid].coord.z;

        umol.atoms[1].charge = molecule.atoms[tid].charge;
        umol.atoms[1].coord.x = molecule.atoms[tid].coord.x;
        umol.atoms[1].coord.y = molecule.atoms[tid].coord.y;
        umol.atoms[1].coord.z = molecule.atoms[tid].coord.z;
        DipoleMoment(umol, u, NULL, NULL, NULL);
        (*h_charge) = molecule.atoms[h_id].charge;
        break;
      }
    }
    else{
      continue;
    }
  }
  DelMolecule(&umol);
}

void FFParamsAssigner(MOLECULE molecule, ForceField *ffield)
{
  size_t i;
  int index;

  for(i = 0; i < molecule.n_atoms; i++){
    index = atomtypeSearch((*ffield), molecule.atoms[i].ainfo.atype_hash, molecule.atoms[i].asymbl, molecule.atoms[i].ainfo.hybrid);
    if(index > -1){
      /* Atom type already in database.
       * Charge recalculation according to the action-value method of reiforcement learning
       * Q(n+1) = Q(n) + 1/n [R(n)-Q(n)] Where:
       * Q(n) is the actual charge of the fingerprint
       * R(n) is the new different charge of the molecule under analysis
       * n is the total number of time that the atomtype was estimated.
       */
      /*
      DEBUG ONLY
      double diff = fabs(ffield->params[index].charge - molecule.atoms[i].charge)/fabs(ffield->params[index].charge);
      printf("%s %s %lu %s %f %f %f\n", molecule.molname, molecule.atoms[i].asymbl, i+1, ffield->params[index].notes, ffield->params[index].charge, molecule.atoms[i].charge, diff);
      */

      ffield->params[index].charge = ffield->params[index].charge + ((1.f/ffield->params[index].freq) * (molecule.atoms[i].charge-ffield->params[index].charge));
      if(strcmp(molecule.atoms[i].asymbl, "H") == 0){
        double u = 0.f;
        double h_charge = 0.f;
        FFHBParamsAssigner(molecule, i, &u, &h_charge);
        ffield->params[index].hb_Umoment = ffield->params[index].hb_Umoment + ((1.f/ffield->params[index].freq) * (u-ffield->params[index].hb_Umoment));
        ffield->params[index].hb_Hcharge = ffield->params[index].hb_Hcharge + ((1.f/ffield->params[index].freq) * (h_charge-ffield->params[index].hb_Hcharge));
      }
      ffield->params[index].freq += 1;

      /* if parameters are != 0
           then
              check if the value are similar or not...
              if not,  do an average... else skip.
        else if parameter are 0 then just write the parameters in index position

      if(ffield->params[index].charge != 0.f){
        double diff = fabs(ffield->params[index].charge - molecule.atoms[i].charge)/fabs(ffield->params[index].charge);
        printf("%s %s %lu %s %f %f %f\n", molecule.molname, molecule.atoms[i].asymbl, i+1, ffield->params[index].notes, ffield->params[index].charge, molecule.atoms[i].charge, diff);
        if(diff > 0.01){
          ffield->params[index].charge = (ffield->params[index].charge + molecule.atoms[i].charge)/2.f;
          ffield->params[index].charge = ffield->params[index].charge + ((1.f/ffield->params[index].freq) * (molecule.atoms[i].charge-ffield->params[index].charge));
          if(strcmp(molecule.atoms[i].asymbl, "H") == 0){
            double u = 0.f;
            double h_charge = 0.f;
            FFHBParamsAssigner(molecule, i, &u, &h_charge);
            ffield->params[index].hb_Umoment += u; ffield->params[index].hb_Umoment /= 2.f;
            ffield->params[index].hb_Hcharge += h_charge; ffield->params[index].hb_Hcharge /= 2.f;
          }
        }
        else{
          continue;
        }
      }
      else{
        ffield->params[index].charge = molecule.atoms[i].charge;
        if(strcmp(molecule.atoms[i].asymbl, "H") == 0){
          double u = 0.f;
          double h_charge = 0.f;
          FFHBParamsAssigner(molecule, i, &u, &h_charge);
          ffield->params[index].hb_Umoment = u;
          ffield->params[index].hb_Hcharge = h_charge;
        }
      }
      */
    }
    else{
      /* Add a new ff parameter */
      ffield->size+=1;
      ffield->params = (Params *) realloc(ffield->params, sizeof(Params)*ffield->size);
      ffield->params[ffield->size-1].atype_hash = strdup(molecule.atoms[i].ainfo.atype_hash);
      ffield->params[ffield->size-1].atom_name = strdup(molecule.atoms[i].asymbl);
      ffield->params[ffield->size-1].hybridtype = molecule.atoms[i].ainfo.hybrid;
      ffield->params[ffield->size-1].charge = molecule.atoms[i].charge;
      ffield->params[ffield->size-1].vdw_R_eq = 0.f; // value to train from experimental data!
      ffield->params[ffield->size-1].vdw_Epsilon = 0.f; // value to train from experimental data!
      ffield->params[ffield->size-1].freq = 1;

      /* if this atom is an hydrogen search for hb_moment and H charge */
      double u = 0.f;
      double h_charge = 0.f;
      if(strcmp(molecule.atoms[i].asymbl, "H") == 0){
        FFHBParamsAssigner(molecule, i, &u, &h_charge);
      }
      else{
        u = h_charge = 0.f;
      }

      ffield->params[ffield->size-1].hb_Umoment = u;
      ffield->params[ffield->size-1].hb_Hcharge = h_charge;
      char buf[21];
      sprintf(buf, "%lu", i+1);
      //ffield->params[ffield->size-1].notes = strdup(molecule.molname);
      ffield->params[ffield->size-1].notes = concat(3, molecule.molname, "_atom_id:", buf);
    }
  }

  /* Check the predicted molecular charge. if not coincident spalm linearly the difference along the atom types
  double inacharge = GetTotalCharge(molecule);
  double ffacharge = 0.f;

  for(i = 0; i < molecule.n_atoms; i++){
    index = atomtypeSearch((*ffield), molecule.atoms[i].ainfo.atype_hash);
    if(index > -1){
       ffacharge += ffield->params[index].charge;
    }
    else{
      printf("Error: atom type %s not found!\n", molecule.atoms[i].ainfo.atype_hash);
    }
  }

  double chdiff = (inacharge - ffacharge)/(double)molecule.n_atoms;

  for(i = 0; i < molecule.n_atoms; i++){
    index = atomtypeSearch((*ffield), molecule.atoms[i].ainfo.atype_hash);
    if(index > -1){
      ffield->params[index].charge += chdiff;
    }
    else{
      printf("Error: atom type %s not found!\n", molecule.atoms[i].ainfo.atype_hash);
    }
  }
  */
}

void SaveFFParams(ForceField ffield, const char *outfile)
{
  size_t i;
  /* Save the forcefield! */
  FILE *f = fopen(outfile, "w");
  if (f == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  /* print the standard text */
  fprintf(f, "# Force field parameters\n");
  fprintf(f, "# 1: atom type hash\n");
  fprintf(f, "# 2: atom name\n");
  fprintf(f, "# 3: hybridtype (0: unspec, 1: s, 2: sp, 3: sp2, 4: sp3, 5: sp3d, 6: sp3d3, 7: ar, 8: other\n");
  fprintf(f, "# 4: atom partial charge\n");
  fprintf(f, "# 5: Van der Waals equilibrium distance (VDW_RStar)\n");
  fprintf(f, "# 6: Van der Waals well depth (VDW_Epsilon)\n");
  fprintf(f, "# 7: Hydrogen bond dipole (HB_dipole)\n");
  fprintf(f, "# 8: Hydrogen bond hydrogen charge (HB_HCharge)\n");
  fprintf(f, "# 9: Atom type frequency\n");
  fprintf(f, "# 10: Atom type notes\n");

  for(i = 0; i < ffield.size; i++){
    //000100000000003001110000000000000000000000000	C 0.0000 0.0000 0.0000 0.0000 0.0000 0) Carbon atom
    fprintf(f, "%s %s %d %.4f %.4f %.4f %.4f %.4f %4lu %s\n",
      ffield.params[i].atype_hash, /* string */
      ffield.params[i].atom_name, /* string */
      ffield.params[i].hybridtype,
      ffield.params[i].charge, /* double */
      ffield.params[i].vdw_R_eq, /* double */
      ffield.params[i].vdw_Epsilon, /* double */
      ffield.params[i].hb_Umoment, /* double */
      ffield.params[i].hb_Hcharge, /* double */
      ffield.params[i].freq, /* size_t */
      ffield.params[i].notes); /* string */
  }
  fclose(f);
}

void DeleteForceField(ForceField *ffield)
{
  for(size_t i = 0; i < ffield->size; i++){
    free(ffield->params[i].atom_name);
    free(ffield->params[i].atype_hash);
    free(ffield->params[i].notes);
  }
  free(ffield->params);
}

/* Binary search algorithm to find the atom type in ForceField
 * Working only with ordered list and not hash list!

int atomtypeSearch(ForceField ffield, char *atype_hash, int begin, int end)
{
  // end is one past the last element, i.e. [begin, end) is a half-open interval.
  if (begin < end){
    int m = (begin+end)/2;
    int cmp = strcmp(ffield.params[m].atype_hash, atype_hash);
    if(cmp == 0)
      return m;
    else if (cmp < 0)
      return atomtypeSearch(ffield, atype_hash, begin, m);
    else
      return atomtypeSearch(ffield, atype_hash, m+1, end);
  }
  else // begin>=end, i.e. no valid array to search
    return -1;
}
*/

int atomtypeSearch(ForceField ffield, char *atype_hash, char *atom_name, int hybridtype)
{
  size_t i;
  int index = -1;
  for(i = 0; i < ffield.size; i++){
    if(strcmp(ffield.params[i].atype_hash, atype_hash) == 0 &&
       strcmp(ffield.params[i].atom_name, atom_name) == 0 &&
       ffield.params[i].hybridtype == hybridtype){
      index = i;
      break;
    }
    else{
      continue;
    }
  }
  return index;
}

void AssignParams(MOLECULE *molecule, ForceField ffield, enum RADIUS_TYPE rtype, int formal_charge)
{
  size_t i;
  int index;
  double dq = 0.f;
  for(i = 0; i < molecule->n_atoms; i++){
    index = atomtypeSearch(ffield, molecule->atoms[i].ainfo.atype_hash, molecule->atoms[i].asymbl, molecule->atoms[i].ainfo.hybrid);
    if(index > -1){
       molecule->atoms[i].ffatompos = index;
       molecule->atoms[i].charge = ffield.params[index].charge;
       dq += molecule->atoms[i].charge;
    }
    else{
      printf("Error: atom type %lu %s %s not found!\n",i+1, molecule->atoms[i].type, molecule->atoms[i].ainfo.atype_hash);
    }

    /*Assign Radius*/
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*Correct the excess of charge*/
  dq -= (double)formal_charge;
  dq /= (double)molecule->n_atoms;
  //printf("dq: %f\n", dq);
  for(i = 0; i < molecule->n_atoms; i++){
    molecule->atoms[i].charge -= dq;
  }
}

void PrintFFParams(ForceField ff)
{
  printf("Atom types :%zu\n", ff.size);
  for(size_t i = 0; i < ff.size; i++){
    printf("%s %s %f %f %f %f %f\n", ff.params[i].atype_hash, ff.params[i].atom_name, ff.params[i].charge, ff.params[i].vdw_R_eq, ff.params[i].vdw_Epsilon, ff.params[i].hb_Umoment, ff.params[i].hb_Hcharge);
  }
}
