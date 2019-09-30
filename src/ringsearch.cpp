/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "ringsearch.h"
#include <algorithm>
#include <vector>
#include <unistd.h>

using namespace std;

RingSearch::RingSearch(MOLECULE molecule)
{
  convert2PathGraph(molecule);
  #ifdef DEBUG
  std::cout << "N. Vertex " << V.size() << " N. Edges " << E.size() << std::endl;
  #endif
  /*
  while(V.size() > 0){
    removeVertex(V[0]);
  }*/

  std::vector <int> V_ = V;
  std::vector< std::vector <int> > E_ = E;
  for(size_t i = 0; i < V_.size(); i++){
    while(V.size() > 0){
      removeVertex(V[0]);
    }
    //shuffle
    for(size_t j = i; j < V_.size(); j++){
      V.push_back(V_[j]);
    }

    for(size_t j = 0; j < i; j++){
      V.push_back(V_[j]);
    }
    E = E_;
  }

  #ifdef DEBUG
  printf("Loops\n");
  for(int i = 0; i < rings.size(); i++){
    for(int j = 0; j < rings[i].size(); j++)
      std::cout << rings[i][j]+1 << "\t";
    std::cout << endl;
  }
  #endif
}

void PrintEdges(std::vector< std::vector<int> > E)
{
  std::cout << "Edge size: " << E.size() << std::endl;
  for(int i = 0; i < E.size(); i++){
    for(int j = 0; j < E[i].size(); j++){
      printf("%d ", E[i][j]+1);
    }
    printf("\n");
  }
  printf("--------------------\n");
}


void RingSearch::branchPruning(std::vector<int> *V, std::vector< std::vector<int> > *E)
{
  std::vector<int> V_;
  for(size_t i = 0; i < (*V).size(); i++){
    size_t c = 0;
    for(size_t j = 0; j < (*E).size(); j++){
      if(std::find((*E)[j].begin(), (*E)[j].end(), (*V)[i]) != (*E)[j].end()){
        c++;
      }
      else{
        continue;
      }
    }
    if(c > 1){
      V_.push_back((*V)[i]);
    }
    else{
      continue;
    }
  }
  size_t prev_E_size = (*E).size();

  //Remove Edges not contained in V_ and replace V with V_
  size_t i = 0;
  while(i < (*E).size()){
    bool rm_edge = false;
    for(size_t j = 0; j < (*E)[i].size(); j++){
      if(std::find(V_.begin(), V_.end(), (*E)[i][j]) != V_.end()){
        continue;
      }
      else{
        rm_edge = true;
        break;
      }
    }

    if(rm_edge == true){
      (*E).erase((*E).begin()+i);
    }
    else{
      i++;
    }
  }

  if(prev_E_size > (*E).size()){
    (*V) = V_;
    branchPruning(V, E);
  }
}

void RingSearch::convert2PathGraph(MOLECULE molecule)
{
  for(size_t i = 0; i < molecule.n_atoms; i++){
    if(molecule.atoms[i].ainfo.connectivity > 1){
      V.push_back(i);
    }
    else{
      continue;
    }
  }

  for(size_t i = 0; i < molecule.n_bonds; i++){
    if(std::find(V.begin(), V.end(), molecule.bonds[i].origin_atom_id) != V.end() &&
       std::find(V.begin(), V.end(), molecule.bonds[i].target_atom_id) != V.end()){
      E.push_back(std::vector<int>());
      E.back().push_back(molecule.bonds[i].origin_atom_id);
      E.back().push_back(molecule.bonds[i].target_atom_id);
    }
    else{
      continue;
    }
  }

  /* Prune vertex with connectivity <= 1 */
  branchPruning(&V,  &E);
}

void join_front_front(std::vector<int> p1, std::vector<int> p2, std::vector<int> *j)
{
  #ifdef DEBUG
  printf("join_front_front a-b a-c = b-a-c\n");
  #endif
  size_t i = p2.size()-1;
  while(i > 0){
    (*j).push_back(p2[i]);
    i--;
  }

  for(int i = 0; i < p1.size(); i++){
    (*j).push_back(p1[i]);
  }
}

void join_front_back(std::vector<int> p1, std::vector<int> p2, std::vector<int> *j)
{
  #ifdef DEBUG
  printf("join_front_back a-b c-a = c-a-b\n");
  #endif
  for(size_t i = 0; i < p2.size()-1; i++){
    (*j).push_back(p2[i]);
  }
  for(size_t i = 0; i < p1.size(); i++){
    (*j).push_back(p1[i]);
  }
}

void join_back_front(std::vector<int> p1, std::vector<int> p2, std::vector<int> *j)
{
  #ifdef DEBUG
  printf("join_back_front b-a a-c = b-a-c\n");
  #endif
  for(size_t i = 0; i < p1.size(); i++){
    (*j).push_back(p1[i]);
  }
  for(size_t i = 1; i < p2.size(); i++){
    (*j).push_back(p2[i]);
  }
}

void join_back_back(std::vector<int> p1, std::vector<int> p2, std::vector<int> *j)
{
  #ifdef DEBUG
  printf("join_back_back b-a c-a = b-a-c\n");
  #endif
  for(size_t i = 0; i < p1.size(); i++){
    (*j).push_back(p1[i]);
  }
  for(size_t i = 0; i < p2.size()-1; i++){
    (*j).push_back(p2[i]);
  }
}

void RingSearch::removeVertex(int id)
{
  #ifdef DEBUG
  printf("Remove %d\n", id+1);
  PrintEdges(E);
  #endif
  int i, j;
  std::vector < std::vector<int> > E_new;
  for(i = 0; i < E.size(); i++){
    if(E[i].front() != E[i].back()){
      for(j = i+1; j < E.size(); j++){
        if(E[j].front() != E[j].back()){
          std::vector<int> join;
          if(E[i].front() == E[j].front() && E[i].front() == id){
            join_front_front(E[i] ,E[j], &join);
          }
          else if(E[i].back() == E[j].front() && E[i].back() == id){
            join_back_front(E[i] ,E[j], &join);
          }
          else if(E[i].front() == E[j].back() && E[i].front() == id){
            join_front_back(E[i] ,E[j], &join);
          }
          else if(E[i].back() == E[j].front() && E[i].back() == id){
            join_back_back(E[i] ,E[j], &join);
          }

          if(join.size() > 0){
            if(isRealPath(join) && isNotPresent(E_new, join) == true){
              E_new.push_back(join);
              #ifdef DEBUG
              printf("Join Case: %d %d %d %d\n", E[i].front()+1, E[i].back()+1, E[j].front()+1, E[j].back()+1);
              PrintEdges(E_new);
              sleep(1);
              #endif
            }
          }
        }
      }
    }
    else{
      continue;
    }
  }

  /* Add the new edges to E */
  for(i = 0; i < E_new.size(); i++){
    E.push_back(E_new[i]);
  }

  i = 0;
  while(i < E.size()){
    if(E[i].front() == E[i].back() && isRealPath(E[i]) && isNotPresent(rings, E[i])){
      rings.push_back(std::vector<int>());
      for(j = 0; j < E[i].size(); j++){
        rings.back().push_back(E[i][j]);
      }
      i++;
    }
    else if(E[i].front() == id || E[i].back() == id){
      E.erase(E.begin()+i);
    }
    else{
      i++;
    }
  }

  #ifdef DEBUG
  std::cout << "Edges for the next loop!" << std::endl;
  PrintEdges(E);
  std::cout << "END" << std::endl;
  #endif

  V.erase(V.begin()+0);
  //V.erase(std::remove(V.begin(), V.end(), id), V.end());
  //PrintEdges(E);
}

bool RingSearch::isRealPath(std::vector<int> path){
  for(size_t i = 1; i < path.size()-1; i++){
    for(size_t j = i+1; j < path.size()-1; j++){
      if(path[i] != path[j]){
        continue;
      }
      else{
        return false;
      }
    }
  }
  return true;
}

bool RingSearch::isNotPresent(std::vector< std::vector<int> > pathgraphs, std::vector<int> path){
  for(size_t i = 0; i < pathgraphs.size(); i++){
    if(compareVectors(pathgraphs[i], path) == true){
      return false;
    }
    else{
      continue;
    }
  }
  return true;
}
