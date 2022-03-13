/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "wfingerprint.h"
#include "molecule.h"
#include "misc.h"
#include "periodic_table.h"

/*
 * Create a string of numbers:
 * 1 = Single bond
 * 2 = Double Bond
 * 3 = Triple Bond
 * 4 = Aromatic
 * 0 = No Bond
 *
 * ATOM   NEAR ATOMS
 * 55100   00000 00000 00000 00000
 * The order is given by the atom weight
 * EXAMPLE... for string of 4
 * [MW1, MW2, MW3, MW4]
 *
 * MW1 > MW2 > MW3 > MW4
 *
 */

typedef struct{
  double mw;
  int type;
  POINT coord;
} ATOMTYPE;

int cmpmw(const void *a, const void *b)
{
  double aa = ((ATOMTYPE *)a)->mw;
  double bb = ((ATOMTYPE *)b)->mw;
  return (aa > bb) ? -1 : (aa < bb) ? 1 : 0;
}


void AtomFragIds(MOLECULE molecule, int aid, int *type, int type_size)
{
//   printf("AID %d\n", aid+1);
  size_t i, j;
  int type_size_ = type_size-1; /*-1 because the last bit in type is the stereoisomerism*/
  ATOMTYPE *atom = malloc(sizeof(ATOMTYPE)*type_size_);
  for(i = 0; i < type_size_; i++){
    atom[i].mw = 0.f;
    atom[i].type = 0;
    atom[i].coord.x = atom[i].coord.y = atom[i].coord.z = 0.f;
    type[i] = 0;
  }
  type[type_size_] = 0;

  char *atomstr;
  j = 0;
  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].target_atom_id == aid){
      if(j < type_size_){
        /*if is a mol2 then split the atom type mol2 to a simple atom*/
        atomstr = strdup(molecule.atoms[molecule.bonds[i].origin_atom_id].type);
        atomstr = strtok(atomstr, ".");
//         printf("Belong to: %d %s %s\n", molecule.bonds[i].origin_atom_id+1, atomstr, molecule.bonds[i].type);
        atom[j].mw = getAtomWeightfromAtomName(atomstr);
        atom[j].coord.x = molecule.atoms[molecule.bonds[i].origin_atom_id].coord.x;
        atom[j].coord.y = molecule.atoms[molecule.bonds[i].origin_atom_id].coord.y;
        atom[j].coord.z = molecule.atoms[molecule.bonds[i].origin_atom_id].coord.z;
        if(strcmp("1", molecule.bonds[i].type) == 0){
          atom[j].type = 1;
        }
        else if(strcmp("2", molecule.bonds[i].type) == 0){
          atom[j].type = 2;
        }
        else if(strcmp("3", molecule.bonds[i].type) == 0){
          atom[j].type = 3;
        }
        else if(strcmp("4", molecule.bonds[i].type) == 0){
          atom[j].type = 4;
        }
        else if(strcmp("ar", molecule.bonds[i].type) == 0){ /*aromatic in MOL2*/
          atom[j].type = 4;
        }
        else{
          atom[j].type = 0;
          printf("Uknown bond type in %s\n", molecule.molname);
        }
//         printf("### %d %f %lu\n", atom[j].type, atom[j].mw, j);
        j++;
        free(atomstr);
      }
      else{
        printf("String out of range in %s", molecule.molname);
      }
    }
    else if(molecule.bonds[i].origin_atom_id == aid){
      if(j < type_size_){
        atomstr = strdup(molecule.atoms[molecule.bonds[i].target_atom_id].type);
        atomstr = strtok(atomstr, ".");
//         printf("Belong to: %d %s %s\n", molecule.bonds[i].target_atom_id+1, atomstr, molecule.bonds[i].type);
        atom[j].mw = getAtomWeightfromAtomName(atomstr);
        atom[j].coord.x = molecule.atoms[molecule.bonds[i].target_atom_id].coord.x;
        atom[j].coord.y = molecule.atoms[molecule.bonds[i].target_atom_id].coord.y;
        atom[j].coord.z = molecule.atoms[molecule.bonds[i].target_atom_id].coord.z;
        if(strcmp("1", molecule.bonds[i].type) == 0){ /*single*/
          atom[j].type = 1;
        }
        else if(strcmp("2", molecule.bonds[i].type) == 0){ /*double*/
          atom[j].type = 2;
        }
        else if(strcmp("3", molecule.bonds[i].type) == 0){ /*triple*/
          atom[j].type = 3;
        }
        else if(strcmp("4", molecule.bonds[i].type) == 0){ /*aromatic in SDF*/
          atom[j].type = 4;
        }
        else if(strcmp("ar", molecule.bonds[i].type) == 0){ /*aromatic in MOL2*/
          atom[j].type = 4;
        }
        /*in sdf there are also 5, 6, 7 and 8 to define bonds.. but here we think that there are only 4 types*/
        else{
          atom[j].type = 0;
          printf("Uknown bond type in %s\n", molecule.molname);
        }
//         printf("### %d %f %lu\n", atom[j].type, atom[j].mw, j);
        j++;
        free(atomstr);
      }
      else{
        printf("String out of range in %s", molecule.molname);
      }
    }
    else{
      continue;
    }
  }

  qsort(atom, type_size_, sizeof(ATOMTYPE), cmpmw);

  /*
  for(i = 0; i < type_size_; i++)
     printf("### %d %f %lu\n", atom[i].type, atom[i].mw, i);

  printf("--------------\n");*/
  /*stereoisomerism function of sdf*/
  type[type_size_] = molecule.atoms[aid].ainfo.stereo;

  /* Check if it's a stereocenter.
   * All the molecular weight must be != 0
   */

  /*
  int stereo = 0;
  int nhyd = 0;
  if(strcmp("C.3", molecule.atoms[aid].type) == 0){
    for(i = 0; i < type_size_; i++){
      if(FLOAT_EQ(atom[i].mw, 0., 1e-6)){
        printf("1 %f\n", atom[i].mw);
        stereo = -1;
        break;
      }
      else if(FLOAT_EQ(atom[i].mw, 1.0079, 1e-6)){
        nhyd++;
      }
      else{
        continue;
      }
    }

    if(stereo == 0)
      if(nhyd > 1)
        stereo = -1;
  }
  else{
    stereo = -1;
  }


  if(stereo == 0){

  */

  /*
    * Determine the Stereocenter R or S
    *
    * Algorithm:
    * http://www.mdpi.org/molecules/papers/61100915/61100915.htm
    *
    * (Xa-Xc) * (Yb-Yc) * (Ya-Yc) * (Xb-Xc)  * Z(H)
    *
    * the sign result correspond to R or S:
    * Negative = R
    * Positive = S
    *
    * Z coord is normally the less heavy atom.. so the last one....
    */

  /*
    double det = (atom[0].coord.x - atom[2].coord.x) *
    (atom[1].coord.y - atom[2].coord.y) *
    (atom[0].coord.y - atom[2].coord.y) *
    (atom[1].coord.x - atom[2].coord.y) * atom[3].coord.z;
    if(det < 0) // R
      type[type_size_] = 5;
    else // S
      type[type_size_] = 6;
  }
  */



  for(i = 0; i < type_size_; i++){
    type[i] = atom[i].type;
//     printf("%f %d\n", atom[i].mw, atom[i].type);
  }
  free(atom);
}

void GetAtomFrag(MOLECULE molecule, int aid, ATOMFRAG *afrag)
{
  size_t i, j;
  /* eacth atom can have 4 bonds and one bit reserved to the stereocenter configuration thus 5 * 5 = 25*/
  int type[5];
  int maxtypesize = 5;
  int maxstrsize = maxtypesize * maxtypesize;

  AtomFragIds(molecule, aid, type, maxtypesize);
  int pos = 0;
  for (j = 0; j < maxtypesize; j++) {
    pos += sprintf(&afrag->str[pos], "%d", type[j]);
  }

  for(i = 0; i < molecule.n_bonds; i++){
    if(molecule.bonds[i].origin_atom_id == aid){
      AtomFragIds(molecule, molecule.bonds[i].target_atom_id, type, maxtypesize);
      for (j = 0; j < maxtypesize; j++) {
        pos += sprintf(&afrag->str[pos], "%d", type[j]);
      }
    }
    else if (molecule.bonds[i].target_atom_id == aid){
      AtomFragIds(molecule, molecule.bonds[i].origin_atom_id, type, maxtypesize);
      for (j = 0; j < maxtypesize; j++) {
        pos += sprintf(&afrag->str[pos], "%d", type[j]);
      }
    }
    else{
      continue;
    }
  }

  for(i = strlen(afrag->str); i < maxstrsize; i++)
    strcat(afrag->str, "0");

  strcpy(afrag->atype, molecule.atoms[aid].type);
}

void ReadWeightTable(char *weightable, ATOMFRAG **afrags, int *afragsize)
{
  FILE *ptr_file;
  char buf[1000];
  char *pch;
  size_t size = 0;
  int pos;

  ptr_file = fopen(weightable,"r");
  if (!ptr_file){
    (*afragsize) = 0;
    return;
  }

  while (fgets(buf,1000, ptr_file) != NULL){
    trim(buf);
    pch = strtok (buf, ",");
    size += 1;
    (*afrags) = realloc((*afrags), sizeof(ATOMFRAG)*size);
    pos = 0;
    while(pch != NULL){
      trim(pch);
      if(pos == 0)
        copystr(pch, (*afrags)[size-1].atype);
      else if(pos == 1)
        copystr(pch, (*afrags)[size-1].str);
      else if(pos == 2)
        (*afrags)[size-1].aweight = atof(pch);
      else
        (*afrags)[size-1].found = atoi(pch);
      pos++;
      pch = strtok (NULL, ",");
    }
  }
  fclose(ptr_file);
  (*afragsize) = size;
}

void InitWeightTable(char *weightable, MOLECULE molecule)
{
  int i, j;
  ATOMFRAG *afrags = NULL;
  ATOMFRAG af;
  FILE *ptr_file;
  int afragsize, ok;
  ReadWeightTable(weightable, &afrags, &afragsize);
  for(i = 0;  i < molecule.n_atoms; i++){
    GetAtomFrag(molecule, i, &af);
    ok = -1;
    for(j = 0; j < afragsize; j++){
      if(strcmp(afrags[j].str, af.str) == 0){
        if(strcmp(afrags[j].atype, af.atype) == 0){
          afrags[j].found += 1;
          ok = 0;
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
    if(ok == 0)
      continue;
    else{
      //Add Fragment to table!
      afragsize += 1;
      afrags = realloc(afrags, sizeof(ATOMFRAG)*afragsize);
      copystr(af.atype, afrags[afragsize-1].atype);
      copystr(af.str, afrags[afragsize-1].str);
      srand((unsigned)afragsize);
      afrags[afragsize-1].aweight = 1;
      afrags[afragsize-1].found = 1;
    }
  }

  ptr_file = fopen(weightable, "w");
  if (!ptr_file)
    return;
  for(i = 0; i < afragsize; i++)
    fprintf(ptr_file,"%s,%s,%f,%d\n", afrags[i].atype, afrags[i].str, afrags[i].aweight, afrags[i].found);
  fclose(ptr_file);
  free(afrags);
}

void CalcPropertFromSumAtomFrag(MOLECULE molecule, char *weightable, double *property)
{
  ATOMFRAG *afrags = NULL;
  ATOMFRAG af;
  int afragsize, i, j, ok;
  ReadWeightTable(weightable, &afrags,  &afragsize);

  (*property) = 0.f;
  for(i = 0;  i < molecule.n_atoms; i++){
    GetAtomFrag(molecule, i, &af);
    ok = -1;
    for(j = 0; j < afragsize; j++){
      if(strcmp(afrags[j].str, af.str) == 0){
        if(strcmp(afrags[j].atype, af.atype) == 0){
          (*property) += afrags[j].aweight;
          ok = 0;
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
    if(ok == 0)
      continue;
    else{
      printf("Fragment not found! \"%s %s\"\n", af.atype, af.str);
    }
  }
  free(afrags);
}

void GetSumAtomFragFP(MOLECULE molecule, char *weightable, double **fp, int *fpsize)
{
  ATOMFRAG *afrags = NULL;
  ATOMFRAG af;
  int afragsize, i, j;
  ReadWeightTable(weightable, &afrags,  &afragsize);
  (*fp) = malloc(sizeof(double)*afragsize);
  for(j = 0; j < afragsize; j++)
    (*fp)[j] = 0.f;
  (*fpsize) = afragsize;
  for(i = 0;  i < molecule.n_atoms; i++){
    GetAtomFrag(molecule, i, &af);
    for(j = 0; j < afragsize; j++){
      if(strcmp(afrags[j].str, af.str) == 0){
        if(strcmp(afrags[j].atype, af.atype) == 0){
          (*fp)[j] += afrags[j].aweight;
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
  }
  free(afrags);
}

void CalcPropertFromProductAtomFrag(MOLECULE molecule, char *weightable, double *property)
{
  ATOMFRAG *afrags = NULL;
  ATOMFRAG af;
  int afragsize, i, j, ok;
  ReadWeightTable(weightable, &afrags,  &afragsize);

  (*property) = 1.;
  for(i = 0;  i < molecule.n_atoms; i++){
    GetAtomFrag(molecule, i, &af);
    ok = -1;
    for(j = 0; j < afragsize; j++){
      if(strcmp(afrags[j].str, af.str) == 0){
        if(strcmp(afrags[j].atype, af.atype) == 0){
          (*property) *= afrags[j].aweight;
          ok = 0;
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
    if(ok == 0)
      continue;
    else{
      printf("Fragment not found! \"%s %s\"\n", af.atype, af.str);
    }
  }
  free(afrags);
}

void WriteAtomFragWeight(MOLECULE molecule, char *weightable, char *outcontribution)
{
  ATOMFRAG *afrags = NULL;
  ATOMFRAG af;
  int afragsize, i, j, ok;
  FILE *ptr_file;
  ptr_file = fopen(outcontribution, "w");
  if(!ptr_file)
    return;

  ReadWeightTable(weightable, &afrags,  &afragsize);

  for(i = 0;  i < molecule.n_atoms; i++){
    GetAtomFrag(molecule, i, &af);
    ok = -1;
    for(j = 0; j < afragsize; j++){
      if(strcmp(afrags[j].str, af.str) == 0){
        if(strcmp(afrags[j].atype, af.atype) == 0){
          fprintf(ptr_file,"%s %s %f\n", afrags[j].atype, afrags[j].str, afrags[j].aweight);
          ok = 0;
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
    if(ok == 0)
      continue;
    else{
      printf("Fragment not found! \"%s %s\"\n", af.atype, af.str);
      fprintf(ptr_file,"%s %s %f\n", afrags[j].atype, afrags[j].str, 0.f);
    }
  }
  free(afrags);
  fclose(ptr_file);
}
