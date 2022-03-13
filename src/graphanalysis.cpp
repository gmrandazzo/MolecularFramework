/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "graphanalysis.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

#ifdef DEBUG
#include <unistd.h>
#endif

int GraphAnalysis::VDistMin(vector<unsigned int> Q, vector<unsigned int> dist){
  unsigned int i, dist_vmin, index_vmin;

  dist_vmin = dist[Q[0]];
  index_vmin = Q[0];

  for(i=0; i<Q.size(); i++){
    if(dist[Q[i]]<dist_vmin){
      dist_vmin = dist[Q[i]];
      index_vmin = Q[i];
    }
  }

  #ifdef DEBUG
  std::cout << "index_v_min " << index_vmin << "\tdist " <<  dist[index_vmin] <<std::endl;
  sleep(1);
  #endif
  return index_vmin;
}

void GraphAnalysis::dijkstra(std::vector< std::vector< int > > graph, int start, vector<int> *node_depth){
  /*
   * 1  function Dijkstra(Graph, source):
   * 2      for each vertex v in Graph:           // Initializations
   * 3          dist[v] := infinity ;              // Unknown distance function from source to v
   * 4          previous[v] := undefined ;         // Previous node in optimal path from source
   * 5      end for ;
   * 6      dist[source] := 0 ;                    // Distance from source to source
   * 7      Q := the set of all nodes in Graph ; // All nodes in the graph are unoptimized - thus are in Q
   * 8      while Q is not empty:                 // The main loop
   * 9          u := vertex in Q with smallest dist[] ;
   * 10          if dist[u] = infinity:
   * 11              break ;                        // all remaining vertices are inaccessible from source
   * 12          fi ;
   * 13          remove u from Q ;
   * 14          for each neighbor v of u:         // where v has not yet been removed from Q.
   * 15              alt := dist[u] + dist_between(u, v) ;
   * 16              if alt < dist[v]:             // Relax (u,v,a)
   * 17                  dist[v] := alt ;
   * 18                  previous[v] := u ;
   * 19              fi  ;
   * 20          end for ;
   * 21      end while ;
   * 22      return dist[] ;
   * 23  end Dijkstra.
   */
  size_t i, j, infinity, u, alt, n_node;

  n_node = graph.size();

  vector<unsigned int> Q(n_node);
  vector<unsigned int> dist(n_node);
  vector<unsigned int> prev(n_node);

  infinity = n_node*n_node;

  for(i = 0; i < n_node; i++){
    Q[i]=i;
    dist[i] = infinity;
    prev[i] = -1;
  }

  dist[start] = 0;

  while(!Q.empty()){

   u = VDistMin(Q, dist);

  if(dist[u] == infinity)
    break;

  //remove vertex u from Q
  for(i = 0; i < Q.size(); i++){
    if(Q[i] == u){
      Q.erase(Q.begin()+i);
    }
  }

  for(i = 0; i < n_node; i++){
    if(graph[u][i] > 0){
      alt = dist[u] + graph[u][i];
      if(alt < dist[i]){
        dist[i] = alt;
        prev[i] = u;
      }
    }
  }
 }

 for(i = 0; i < n_node; i++){
   (*node_depth).push_back(dist[i]);
 }
}
