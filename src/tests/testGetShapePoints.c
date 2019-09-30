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
#include "../atomdesc.h"
#include "../miscdesc.h"
#include "../shapedesc.h"
#include "../geomdesc.h"

int main(int argc, char **argv)
{
  if(argc != 3){
    printf("\nUsage %s file.mol2 outshape.mol2\n\n", argv[0]);
    return -1;
  }

  SHAPEPNT *pnt;
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  getShapePoints(molecule, &pnt, vanderwaals);
  WriteMol2SurfPoints(pnt, argv[2]);
  PrintMolecule(molecule);
  deleteShapePoints(&pnt);
  DelMolecule(&molecule);
  return 0;
}
