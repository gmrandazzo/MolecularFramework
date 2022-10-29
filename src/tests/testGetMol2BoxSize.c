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

  double cx, cy, cz, xmin, xmax, ymin, ymax, zmin, zmax;
  xmin = xmax = molecule.atoms[0].coord.x;
  ymin = ymax = molecule.atoms[0].coord.y;
  zmin = zmax = molecule.atoms[0].coord.z;
  cx = cy = cz = 0.;
  for(int i = 0; i < molecule.n_atoms; i++){
    cx += molecule.atoms[i].coord.x;
    cy += molecule.atoms[i].coord.y;
    cz += molecule.atoms[i].coord.z;
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
  
  cx/=(float)molecule.n_atoms;
  cy/=(float)molecule.n_atoms;
  cz/=(float)molecule.n_atoms;
  printf("xmin ymin zmin: %f %f %f\n", xmin, ymin, zmin);
  printf("xmax ymax zmax: %f %f %f\n", xmax, ymax, zmax);
  printf("xsize ysize zsize: %f %f %f\n", xmax-xmin, ymax-ymin, zmax-zmin);
  printf("center: %f %f %f\n", cx, cy, cz);

  //GetMolBox(, double size, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax)

  DelMolecule(&molecule);
  return 0;
}
