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
#include "../massanalysis.h"
#include "../shapedesc.h"
#include "../geomdesc.h"

int main(int argc, char **argv)
{
  if(argc == 1){
    printf("\nUsage %s file.mol2\n\n", argv[0]);
    return -1;
  }

  char *bform;
  double planarity, /*lenght,*/ bfmw, exbfmw, dm, volume, surface;
  int nm, n_non_H, n_C, n_H, n_N, n_P, n_S, n_O, n_F, n_Cl, n_Br, n_I;
  MOLECULE molecule;
  NewMOL2Molecule(&molecule, argv[1]);
  AtomAnalyzer(&molecule, 1);
  /* Print tests
  PrintMolecule(molecule);
  PrintMoleculeAtomInfo(molecule);
  PrintMoleculeAtomFingerprint(molecule);*/

  /* Descriptor tests */
  GetMolecularFormula(molecule, &bform);
  GetPlanarity(molecule, &planarity);
  /*GetMoleculaLenght(molecule, &lenght);*/
  GetMWfromMolecularFormula(bform, &bfmw);
  GetExactMWfromMolecularFormula(bform, &exbfmw);
  GetNominalMW(molecule, &nm);
  DipoleMoment(molecule, &dm, NULL, NULL, NULL);
  /*GetMW(molecule, &bfmw);*/
  GetVDWMolVolSurf(molecule, &volume, &surface, 2.0);
  n_non_H = GetNNonHydrogenAtoms(molecule);
  n_C = GetNCarbonAtoms(molecule);
  n_H = GetNHydrogenAtoms(molecule);

  n_N =GetNNitrogenAtoms(molecule);
  n_O = GetNOxygenAtoms(molecule);
  n_S = GetNSulfurAtoms(molecule);
  n_P = GetNPhosphorousAtoms(molecule);
  n_F = GetNFluorineAtoms(molecule);
  n_Cl = GetNChlorineAtoms(molecule);
  n_Br = GetNBromineAtoms(molecule);
  n_I = GetNIodineAtoms(molecule);

  printf("%s %s %.3f %.10f %.10f %d %.3f %d %d %d %d %d %d %d %d %d %d %d %.3f %.3f\n",
  molecule.molname, bform, planarity, bfmw, exbfmw, nm, dm,
  n_non_H, n_C, n_H, n_N, n_O, n_S, n_P, n_F, n_Cl, n_Br, n_I, volume, surface);

  DelMolecule(&molecule);
  free(bform);
  return 0;
}
