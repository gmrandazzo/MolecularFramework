/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../atomanalysis.h"
#include "../molecule.h"
#include "../cppwrap.h"
#include "../miscdesc.h"


int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2\n\n", argv[0]);
    return -1;
  }
  printf("%s,", argv[1]);
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);

  double xmin, xmax, ymin, ymax, zmin, zmax;
  xmin = xmax = molecule.atoms[0].coord.x;
  ymin = ymax = molecule.atoms[0].coord.y;
  zmin = zmax = molecule.atoms[0].coord.z;
  for(int i = 0; i < molecule.n_atoms; i++){
    if(molecule.atoms[i].coord.x < xmin)
      xmin = molecule.atoms[i].coord.x;
    if(molecule.atoms[i].coord.x > xmax)
      xmax = molecule.atoms[i].coord.x;
    if(molecule.atoms[i].coord.y < ymin)
      ymin = molecule.atoms[i].coord.y;
    if(molecule.atoms[i].coord.y > ymax)
      ymax = molecule.atoms[i].coord.y;
    if(molecule.atoms[i].coord.z < zmin)
      zmin = molecule.atoms[i].coord.z;
    if(molecule.atoms[i].coord.z > zmax)
      zmax = molecule.atoms[i].coord.z;
  }

  printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n", xmin, xmax, ymin, ymax, zmin, zmax, xmax-xmin, ymax-ymin, zmax-zmin);

  //GetMolBox(, double size, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax)

  DelMolecule(&molecule);
  return 0;
}
