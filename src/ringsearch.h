/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef RINGSEARCH_H
#define RINGSEARCH_H
#include <iostream>
#include <vector>

using namespace std;

#include "molecule.h"

#include <algorithm>

/*
 * A New Algorithm for Exhaustive Ring Perception in a Molecular Graph
 * J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152
 * DOI: 10.1021/ci960322f
 */
class RingSearch
{
public:
  RingSearch(MOLECULE molecule);
  std::vector< std::vector<int> > getRings(){ return rings; }

private:
  std::vector< std::vector<int> > rings;
  std::vector<int> V;
  std::vector< std::vector<int> > E;
  void convert2PathGraph(MOLECULE molecule);
  void branchPruning(std::vector<int> *V, std::vector< std::vector<int> > *E);
  void removeVertex(int id);
  bool isRealPath(std::vector<int> path);
  bool isNotPresent(std::vector< std::vector<int> > pathgraphs, std::vector<int> path);
  template <class T>
  static bool compareVectors(std::vector<T> a, std::vector<T> b)
  {
    if (a.size() != b.size()){
      return false;
    }
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return (a == b);
  }
};

#endif
