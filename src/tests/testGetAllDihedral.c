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
    double dhangle=0.f;
    uivector *dhids;
    NewMOL2Molecule(&m, argv[1]);

    for(int i = 0; i < m.n_atoms; i++){
      for(int j = i+1; j < m.n_atoms; j++){
        initUIVector(&dhids);
        dhangle = 0.f;
        Get3DDihedralAngle(m,
                           j,
                           i,
                           &dhangle,
                           &dhids);
       if(FLOAT_EQ(dhangle, -9999.f, 1e-1) != 1){
         printf("ids: %zu %zu %zu %zu angle: %f\n", dhids->data[0],
                                                    dhids->data[1],
                                                    dhids->data[2],
                                                    dhids->data[3],
                                                    fabs(Rad2Grad(dhangle)));
       }
       DelUIVector(&dhids);
      }
    }
    DelMolecule(&m);

  }
  return 0;
}
