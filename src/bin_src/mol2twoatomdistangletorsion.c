/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../molecule.h"
#include "../geomdesc.h"

int main(int argc, char **argv)
{
  if(argc != 4){
    printf("\nUsage %s [file.mol2] [atom_id1] [atom_id2]\n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);
    int atom_id1 = atoi(argv[2]);
    int atom_id2 = atoi(argv[3]);
    double dst = GetDistance_(molecule.atoms[atom_id1].coord, molecule.atoms[atom_id2].coord);
    double angle, dihedral;
    Get3DAngle(molecule, atom_id1, atom_id2, &angle, NULL);
    Get3DDihedralAngle(molecule, atom_id1, atom_id2, &dihedral, NULL);

    if(FLOAT_EQ(angle, -9999.0, 1e-3)){
      /*No angle possible*/
      angle = -9999.0;
    }
    else{
      angle = fabs(Rad2Grad(angle));
    }

    if(FLOAT_EQ(dihedral, -9999.0, 1e-3)){
      /* No dihedral possible */
      dihedral = -9999.0;
    }
    else{
      dihedral = fabs(Rad2Grad(dihedral));
    }

    printf("%s_%d%s->%d%s,%f,%f,%f\n", molecule.molname,
                                       atom_id1,
                                       molecule.atoms[atom_id1].type,
                                       atom_id2,
                                       molecule.atoms[atom_id2].type,
                                       dst,
                                       angle,
                                       dihedral);
    DelMolecule(&molecule);
    return 0;
  }
}
