/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef GRAPHANALYSIS_H
#define GRAPHANALYSIS_H
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class GraphAnalysis
{
public:
  GraphAnalysis(){}
  void dijkstra(std::vector< std::vector< int > > graph, int start, vector<int> *node_depth);

private:
  /* Utilised by dijkstra */
  int VDistMin(vector<unsigned int> Q, vector<unsigned int> dist);

};

#endif
