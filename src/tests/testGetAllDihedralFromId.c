/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../molecule.h"
#include "../geomdesc.h"

int main(int argc, char **argv)
{
  if(argc != 3){
    printf("\nUsage %s [mol2] [atom id]\n\n", argv[0]);
  }
  else{
    MOLECULE m;
    double dhangle=0.f;
    uivector *dhids;
    NewMOL2Molecule(&m, argv[1]);
    int atom_id = atoi(argv[2]);



    printf("Dihedral starging from %d\n", atom_id);
    for(int i = 0; i < m.n_atoms; i++){
      if(i == atom_id){
        continue;
      }
      else{
        initUIVector(&dhids);
        Get3DDihedralAngle(m,
                           atom_id,
                           i,
                           &dhangle,
                           &dhids);
       if(FLOAT_EQ(dhangle, -9999.f, 1e-1) != 1){
         printf("ids: %zu %zu %zu %zu angle: %f\n", dhids->data[0],
                                                    dhids->data[1],
                                                    dhids->data[2],
                                                    dhids->data[3],
                                                    dhangle);

       }
       DelUIVector(&dhids);
      }
    }
    DelMolecule(&m);

  }
  return 0;
}
