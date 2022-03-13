/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "kekule.h"

#include <iostream>
#include <vector>


typedef struct
{
  std::vector<int> neighbours;
  size_t valence;
  size_t lonepairs;
  size_t aid;
} KATOM;


void Kekulize(MOLECULE *molecule);
