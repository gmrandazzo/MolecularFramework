/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <math.h>
#include "../molecule.h"
#include "../geomdesc.h"

int main(int argc, char **argv)
{
  if(argc != 2){
    printf("\nUsage %s [mol2]\n\n", argv[0]);
  }
  else{
    MOLECULE m;
    double angle=0.f;
    uivector *aids;
    NewMOL2Molecule(&m, argv[1]);

    for(int i = 0; i < m.n_atoms; i++){
      for(int j = i+1; j < m.n_atoms; j++){
        initUIVector(&aids);
        Get3DAngle(m,
                   j,
                   i,
                   &angle,
                   aids);
       if(FLOAT_EQ(angle, -9999.f, 1e-1) != 1){
         printf("ids: %zu %zu %zu angle: %f\n", aids->data[0],
                                                aids->data[1],
                                                aids->data[2],
                                                fabs(Rad2Grad(angle)));

       }
       DelUIVector(&aids);
      }
    }
    DelMolecule(&m);

  }
  return 0;
}
