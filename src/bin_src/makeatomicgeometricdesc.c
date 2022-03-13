/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include "../molecule.h"
#include "../geomdesc.h"


int main(int argc, char **argv)
{
  if(argc != 2){
    printf("\nUsage %s [file.mol2] \n\n", argv[0]);
    return -1;
  }
  else{
    int i, j;
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    
    for(i = 0; i < molecule.n_atoms; i++){
      for(j = 0; j < molecule.n_bond; j++){
        /* Calculates the angle (in radians) between two vectors pointing outward from one center */
        double GetDistance_(POINT p1, POINT p2);

        /*Calculate angle between p1-c-p2 in radians */
        double GetAngle(POINT p1, POINT p2, POINT c);
      }
    }

    DelMolecule(&molecule);
    return 0;
  }
}
