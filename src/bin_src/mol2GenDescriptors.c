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
#include "../massanalysis.h"
#include "../miscdesc.h"
#include "../shapedesc.h"


int main(int argc, char **argv)
{
  if(argc != 2){
    printf("\nUsage %s [file.mol2] \n\n", argv[0]);
    return -1;
  }
  else{
    MOLECULE molecule;
    NewMOL2Molecule(&molecule, argv[1]);


    double mw, exmw, volume, surface, diffusion, planarity;

    GetPlanarity(molecule, &planarity);
    GetMW(molecule, &mw);
    GetExactMW(molecule, &exmw);
    GetVDWMolVolSurf(molecule, &volume, &surface, 0.5);
    CalcMolecularDiffusion_StockesEinstein(molecule, &diffusion);

    printf("%s,", molecule.molname);
    printf("%f,", planarity);
    printf("%f,", mw);
    printf("%f,", exmw);
    printf("%f,", volume);
    printf("%f,", surface);
    printf("%f,", diffusion);
    printf("%f,", ElectronicChemicalPotentialFromMolecule(molecule));
    printf("%f,", GetSumFirstIonizationPotetial(molecule));
    printf("%d,", GetRotableBonds(molecule));
    printf("%d,", GetNCarbonAtoms(molecule));
    printf("%d,", GetNSP3CarbonAtoms(molecule));
    printf("%d,", GetNSP2CarbonAtoms(molecule));
    printf("%d,", GetNSPCarbonAtoms(molecule));
    printf("%d,", GetNArCarbonAtoms(molecule));
    printf("%d,", GetNCatCarbonAtoms(molecule));
    printf("%d,", GetNHydrogenAtoms(molecule));
    printf("%d,", GetNNitrogenAtoms(molecule));
    printf("%d,", GetNSP3NitrogenAtoms(molecule));
    printf("%d,", GetNPlanarNitrogenAtoms(molecule));
    printf("%d,", GetNSP2NitrogenAtoms(molecule));
    printf("%d,", GetNSPNitrogenAtoms(molecule));
    printf("%d,", GetNArNitrogenAtoms(molecule));
    printf("%d,", GetNamNitrogenAtoms(molecule));
    printf("%d,", GetNOxygenAtoms(molecule));
    printf("%d,", GetNSP3OxygenAtoms(molecule));
    printf("%d,", GetNSP2OxygenAtoms(molecule));
    printf("%d,", GetNCO2OxygenAtoms(molecule));
    printf("%d,", GetNSulfurAtoms(molecule));
    printf("%d,", GetNSP3SulfurAtoms(molecule));
    printf("%d,", GetSOSulfurAtoms(molecule));
    printf("%d,", GetNSO2SulfurAtoms(molecule));
    printf("%d,", GetNPhosphorousAtoms(molecule));
    printf("%d,", GetNFluorineAtoms(molecule));
    printf("%d,", GetNChlorineAtoms(molecule));
    printf("%d,", GetNBromineAtoms(molecule));
    printf("%d,", GetNIodineAtoms(molecule));
    printf("%d,", GetNNonHydrogenAtoms(molecule));
    printf("%d,", GetNOHBonds(molecule));
    printf("%d,", GetNNHBonds(molecule));
    printf("%d,", GetNSHBonds(molecule));
    printf("%d,", GetNPHBonds(molecule));
    printf("%d,", GetCCDoubleBonds(molecule));
    printf("%d,", GetCODoubleBonds(molecule));
    printf("%d,", GetNNNDoubleBonds(molecule));
    printf("%d,", GetNNODoubleBonds(molecule));
    printf("%d,", GetNSODoubleBonds(molecule));
    printf("%d,", GetNPODoubleBonds(molecule));
    printf("%d,", GetCCTripleBonds(molecule));
    printf("%d,", GetNCNTripleBonds(molecule));
    printf("%d\n", GetSumValenceElectrons(molecule));
    DelMolecule(&molecule);
    return 0;
  }
}
