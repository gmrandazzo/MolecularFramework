/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../molecule.h"
#include "../massanalysis.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s <BRUTE FORMULA> <0 for mass exact | 1 for molecular weight>\n\n", argv[0]);
    return -1;
  }

  double mass;

  if(atoi(argv[2]) == 0){
    GetExactMWfromMolecularFormula(argv[1], &mass);
  }
  else{
    GetMWfromMolecularFormula(argv[1], &mass);
  }

  printf("%s %.10f\n", argv[1], mass);
  return 0;
}
