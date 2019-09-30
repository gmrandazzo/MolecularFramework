/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../atomanalysis.h"
#include "../geomdesc.h"
#include "../molecule.h"

/*
void ReadCoord(char *f, matrix **coord)
{
  FILE *ptr_file;
  dvector *row;
  char buf[1000];
  NewDVector(&row, 3);

  ptr_file =fopen(f, "r");
  if(!ptr_file)
    return;

  while(fgets(buf, sizeof(buf), ptr_file) != NULL){
    sscanf(buf, "%lf %lf %lf", &row->data[0], &row->data[1], &row->data[2]);
    MatrixAppendRow(coord, row);
  }

  fclose(ptr_file);
  DelDVector(&row);
}
*/

int main(int argc, char **argv)
{
  if(argc > 1){
    int i;
    double planarity;
    for(i = 1; i < argc; i++){
      MOLECULE molecule;
      NewMOL2Molecule(&molecule, argv[i]);
      AtomAnalyzer(&molecule, 1);
      GetPlanarity(molecule, &planarity);
      printf("%s\t%lf\n", molecule.molname, planarity);
      DelMolecule(&molecule);
    }
  }
  else{
    printf ("\nUsage: %s file1.mol2 file2.mol2 file3.mol2 ... fileN.mol2 \n\n", argv[0]);
  }
  return 0;
}
