/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */


/*#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <ForceField/ForceField.h>
#include <ForceField/UFF/Params.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
*/

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>

//#include "RingSearch.h"
#include "ringsearch.h"

#include "cppwrap.h"
#include "molecule.h"
#include "graphanalysis.h"
#include "geomanalysis.h"
#include "scientific.h"


using namespace std;

void rdkitMol2Read(const char *filename, MOLECULE *molecule)
{
  // Get other information using RDKit
  // MMFF R_ii*, epsilon_ii for MIFs

  /*RDKit::RWMol *mol = (RDKit::RWMol*)0;
  try {
    rdErrorLog->df_enabled = false; // Do not report parsing errors
    mol = RDKit::Mol2FileToMol(filename, 1, 0);
    //ForceFields::UFF::UFFVdW uffVdWParams;
    for(size_t i = 0; i < mol->getNumAtoms(); i++){


      //Atom.h
      //  typedef enum {
      //  81     UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
      //  82     S,
      //  83     SP,
      //  84     SP2,
      //  85     SP3,
      //  86     SP3D,
      //  87     SP3D2,
      //  88     OTHER  //!< unrecognized hybridization
      //  89   } HybridizationType;
      //  90

      molecule->atoms[i].ainfo.hybrid = (hybridtype) mol->getAtomWithIdx(i)->getHybridization();

      // Get aromaticity
      if(mol->getAtomWithIdx(i)->getIsAromatic() == true){
        molecule->atoms[i].ainfo.aromatic = 1;
      }
      else{
        molecule->atoms[i].ainfo.aromatic = 0;
      }

      // Get stereochemistry

      //Atom.h
      //    typedef enum {
      //93     CHI_UNSPECIFIED = 0,  //!< chirality that hasn't been specified
      //94     CHI_TETRAHEDRAL_CW,   //!< tetrahedral: clockwise rotation (SMILES \@\@)
      //95     CHI_TETRAHEDRAL_CCW,  //!< tetrahedral: counter-clockwise rotation (SMILES
      //96                           //\@)
      //97     CHI_OTHER             //!< some unrecognized type of chirality
      //98   } ChiralType;

      molecule->atoms[i].ainfo.stereo = mol->getAtomWithIdx(i)->getChiralTag();


      //RDKit::UFF::getUFFVdWParams(*mol, i, i, uffVdWParams);
      //molecule->atoms[i].vdw_Rstar_i = uffVdWParams.x_ij;
      //molecule->atoms[i].vdw_Eps_i = uffVdWParams.D_ij;
    }
    delete mol;
  }
  catch (...) { // RDKit::MolSanitizeException and friends
    mol = (RDKit::RWMol*)0;
    delete mol;
  }*/
}

void ringSearch(MOLECULE *molecule)
{
  //from ringsearch.cpp
  for(size_t i = 0; i < (*molecule).n_atoms; i++)
  (*molecule).atoms[i].ainfo.cycle = 0;

  RingSearch rs((*molecule));
  std::vector< std::vector<int> > rings =  rs.getRings();

  if((*molecule).rings != NULL){
   for(size_t i = 0; i < (*molecule).n_rings; i++){
     free((*molecule).rings[i].atoms);
   }
   free((*molecule).rings);
  }
  (*molecule).rings = (RINGS*)malloc(sizeof(RINGS)*rings.size());

  (*molecule).n_rings = rings.size();
  for(size_t i = 0; i < rings.size(); i++){
    (*molecule).rings[i].atoms = (int*)malloc(sizeof(int)*rings[i].size());
    /*avoid to add the last item since is an id repetition*/
    (*molecule).rings[i].size = rings[i].size()-1;
    for(size_t j = 0; j < rings[i].size()-1; j++){
      (*molecule).atoms[rings[i][j]].ainfo.cycle = 1;
      (*molecule).rings[i].atoms[j] = rings[i][j];
    }
  }
}

/* Create fingerprint for atomtype definition */
 enum position{
   c_sp = 0,
   c_sp2,
   c_sp3,
   c_ar,
   n_sp,
   n_sp2,
   n_sp3,
   n_ar,
   o_sp2,
   o_sp3,
   o_ar,
   s_sp2,
   s_sp3,
   s_sp3d,
   s_sp3d2,
   s_ar,
   p,
   f,
   cl,
   br,
   i,
   h,
   other
 };

 struct enumtypes
 {
    string str;
    position pos;
 };

 struct enumtypes atypes[] = {
   {"C", c_sp},
   {"C", c_sp2},
   {"C", c_sp3},
   {"C", c_ar},
   {"N", n_sp},
   {"N", n_sp2},
   {"N", n_sp3},
   {"N", n_ar},
   {"O", o_sp2},
   {"O", o_sp3},
   {"O", o_ar},
   {"S", s_sp2},
   {"S", s_sp3},
   {"S", s_sp3d},
   {"S", s_sp3d2},
   {"S", s_ar},
   {"P", p},
   {"F", f},
   {"Cl", cl},
   {"Br", br},
   {"I", i},
   {"H", h},
   {"Other", other} // 23
 };

static int GetFPPosition(MOLECULE *molecule, int k){
  char *atomsym  = molecule->atoms[k].asymbl;
  //strcpy(atomsym, molecule->atoms[k].asymbl);
  if(strcmp(atomsym, "C") == 0 && molecule->atoms[k].ainfo.hybrid == SP){
    return c_sp;
  }
  else if(strcmp(atomsym, "C") == 0 && molecule->atoms[k].ainfo.hybrid == SP2){
      return c_sp2;
  }
  else if(strcmp(atomsym, "C") == 0 && molecule->atoms[k].ainfo.hybrid == SP3){
    /* 0: Nochiral 1: R, 2: S, 3: ar 4: OTHER? */
    //if(molecule->atoms[k].ainfo.stereo == 0)
    return c_sp3;
  }
  else if(strcmp(atomsym, "C") == 0 && molecule->atoms[k].ainfo.hybrid == AR){
      return c_ar;
  }
  else if(strcmp(atomsym, "N") == 0 && molecule->atoms[k].ainfo.hybrid == SP){
    return n_sp;
  }
  else if(strcmp(atomsym, "N") == 0 && molecule->atoms[k].ainfo.hybrid == SP2){
    return n_sp2;
  }
  else if(strcmp(atomsym, "N") == 0 && molecule->atoms[k].ainfo.hybrid == SP3){
    return n_sp3;
  }
  else if(strcmp(atomsym, "N") == 0 && molecule->atoms[k].ainfo.hybrid == AR){
    return n_ar;
  }
  else if(strcmp(atomsym, "O") == 0 && molecule->atoms[k].ainfo.hybrid == SP2){
    return o_sp2;
  }
  else if(strcmp(atomsym, "O") == 0 && molecule->atoms[k].ainfo.hybrid == SP3){
    return o_sp3;
  }
  else if(strcmp(atomsym, "O") == 0 && molecule->atoms[k].ainfo.hybrid == SP3){
    return o_sp3;
  }
  else if(strcmp(atomsym, "P") == 0){
    return p;
  }
  else if(strcmp(atomsym, "S") == 0 && molecule->atoms[k].ainfo.hybrid == SP2){
    return s_sp2;
  }
  else if(strcmp(atomsym, "S") == 0 && molecule->atoms[k].ainfo.hybrid == SP3){
    return s_sp3;
  }
  else if(strcmp(atomsym, "S") == 0 && molecule->atoms[k].ainfo.hybrid == SP3D){
    return s_sp3d;
  }
  else if(strcmp(atomsym, "S") == 0 && molecule->atoms[k].ainfo.hybrid == SP3D2){
    return s_sp3d2;
  }
  else if(strcmp(atomsym, "F") == 0){
    return f;
  }
  else if(strcmp(atomsym, "Cl") == 0){
    return cl;
  }
  else if(strcmp(atomsym, "Br") == 0){
    return br;
  }
  else if(strcmp(atomsym, "I") == 0){
    return i;
  }
  else if(strcmp(atomsym, "H") == 0){
    return h;
  }
  else{
    return other;
  }
  //free(atomsym);
}

template <typename T>
  std::string NumberToString ( T Number )
  {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }

void AdjMatGen(MOLECULE *molecule, std::vector< std::vector<int> > *adj_mx){
 int i, j;
 for(i = 0; i < molecule->n_atoms; i++){
   (*adj_mx)[i][i] = 0;
   for(j = 0; j < molecule->n_bonds; j++){
     if(molecule->bonds[j].origin_atom_id == i){
       (*adj_mx)[i][molecule->bonds[j].target_atom_id] = 1;
     }
     else if(molecule->bonds[j].target_atom_id == i){
       (*adj_mx)[i][molecule->bonds[j].origin_atom_id] = 1;
     }
     else
      continue;
   }
 }
}


struct VEDGE {
  VEDGE(int o_, int t_)
  {
    o = o_;
    t = t_;
  }
  int o, t;
};

int edgeIsVisited(std::vector<VEDGE> visited, int o, int t)
{
  int i;
  for(i = 0; i < visited.size(); i++){
    if(visited[i].o == o && visited[i].t == t){
      return 0;
    }
    else{
      continue;
    }
  }
  return -1;
}


void getAtomDepth(MOLECULE *molecule, int src, int depth_max, int **atom_depth_)
{
  int i, j, k, depth_max_;
  /*  Calculate atom depth from the i atoms randazzo version*/
  std::vector < VEDGE > v;
  std::vector < int > prev_idx;

  std::vector<int> atom_depth;
  for(j = 0; j < molecule->n_atoms; j++){
    atom_depth.push_back(9999);
  }

  /*if(depth_max > molecule->n_atoms || depth_max < 0){
    depth_max_ = molecule->n_atoms;
  }*/
  depth_max_ = depth_max;

  atom_depth[src] = 0;
  int depth = 1; // We start from depth 1
  prev_idx.push_back(src);
  while(depth < depth_max_/*depth_max+1*/){

    //std::cout << depth << std::endl;

    std::vector < int > next_idx;
    for(k = 0; k < prev_idx.size(); k++){
      for(j = 0; j < molecule->n_bonds; j++){
        if(molecule->bonds[j].origin_atom_id == prev_idx[k] &&
          edgeIsVisited(v, molecule->bonds[j].origin_atom_id, molecule->bonds[j].target_atom_id) == -1) {
          next_idx.push_back(molecule->bonds[j].target_atom_id);
          atom_depth[molecule->bonds[j].target_atom_id] = depth;
          v.push_back(VEDGE(molecule->bonds[j].origin_atom_id, molecule->bonds[j].target_atom_id));
        }
        else if(molecule->bonds[j].target_atom_id == prev_idx[k] &&
          edgeIsVisited(v, molecule->bonds[j].origin_atom_id, molecule->bonds[j].target_atom_id) == -1){
          next_idx.push_back(molecule->bonds[j].origin_atom_id);
          atom_depth[molecule->bonds[j].origin_atom_id] = depth;
          v.push_back(VEDGE(molecule->bonds[j].origin_atom_id, molecule->bonds[j].target_atom_id));
        }
        else{
          continue;
        }
      }
    }
    //prev_idx = next_idx;
    prev_idx.assign(next_idx.begin(), next_idx.end());
    depth++;
  }

  (*atom_depth_) = (int*) malloc(sizeof(int)*atom_depth.size());
  for(i = 0; i < atom_depth.size(); i++){
    (*atom_depth_)[i] = atom_depth[i];
  }
  //return atom_depth.data();
}

void atomToFingerPrint(MOLECULE *molecule, int depth_max)
{
  int i, j, k, l;
  /*
  // allocate adjacence matrix for dijkstra version
  std::vector< vector<int> > adj_mx;

  for(i = 0; i < molecule->n_atoms; i++){
    adj_mx.push_back(std::vector<int>());
    for(j = 0; j < molecule->n_atoms; j++)
      adj_mx[i].push_back(0);
  }

  // Create adjacence matrix with same weight for dijkstra
  AdjMatGen(molecule, &adj_mx);
  */

  int fpsize = sizeof(atypes)/sizeof(atypes[0]);
  GraphAnalysis ga;

  for(i = 0; i < molecule->n_atoms; i++){
    /* randazzo version*/
    int *atom_depth;
    getAtomDepth(molecule, i, depth_max+1, &atom_depth);

    /* dijkstra version
    vector<int> atom_depth;
    ga.dijkstra(adj_mx, i, &atom_depth);
    */

    /* Create bitstring */
    std::vector< std::vector<int> > fp_str;
    for(j = 0; j < depth_max; j++){
      fp_str.push_back(std::vector<int>());
      for(k = 0; k < fpsize; k++){
        fp_str[j].push_back(0);
      }
    }

    for(j = 1; j < depth_max+1; j++){ /* Start from 1 bond lenght and not from 0 */
      for(k = 0; k < molecule->n_atoms; k++){
        if(atom_depth[k] == j){
          int pos = GetFPPosition(molecule, k);
          if(pos != -1){
            fp_str[j-1][pos] += 1;
          }
          else{
            fp_str[j-1][fpsize-1] += 1;
          }
        }
        else{
          continue;
        }
      }
    }
    // Create the atom fingerprint string
    std::string sum_fpstr;
    for(j = 0; j < depth_max; j++){
      for(k = 0; k < fpsize; k++){
        if(j == 0 && k == 0){
          sum_fpstr.append(NumberToString(fp_str[j][k]));
        }
        else{
          sum_fpstr.append(","+NumberToString(fp_str[j][k]));
        }
      }
    }
    //std::cout << "Size string: "<< sum_fpstr.size() << std::endl;
    free(atom_depth);
    // Copy the fingerprint output
    molecule->atoms[i].ainfo.atype_hash = strdup(sum_fpstr.c_str());
  }
}


/*void Align3DConformations(MOLECULE *from, MOLECULE *to)
{
  GeomAnalysis ga;
  ga.SimplexAlign3DConformations(from, to);
}
*/
