/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "periodic_table.h"

typedef struct{
   char *atom;
   double atomval;
} PeriodicTable;

/*
 * Periodic Table according to
 * https://www.degruyter.com/downloadpdf/j/pac.2016.88.issue-3/pac-2015-0305/pac-2015-0305.pdf
 * some elements verified with http://www.chemcalc.org/
 * and organized as:
 * Atom name, Relative atomic mass
 */
PeriodicTable relamptable[] = {
  {"H", 1.007941},
  {"D", 2.01410178}, /*Deuterium isotope */
  {"Li", 	6.940038},
  {"Na", 22.989769282},
  {"K", 39.098301},
  {"Rb", 85.467664},
  {"Cs", 132.90545192},
  {"Fr", 223.02},
  {"Be", 9.0121823},
  {"Mg", 24.305052},
  {"Ca", 40.0784},
  {"Sr", 87.621},
  {"Ba", 137.3277},
  {"Ra", 226},
  {"B", 10.8117},
  {"Al", 26.98153868},
  {"Ga", 69.7231},
  {"In", 114.8183},
  {"Tl", 204.38332},
  {"C", 12.010736},
  {"Si", 28.085499},
  {"Ge", 72.631},
  {"Sn", 118.7107},
  {"Pb", 207.21},
  {"N", 14.006703},
  {"P", 30.9737622},
  {"As", 74.921602},
  {"Sb", 121.7601},
  {"Bi", 208.980401},
  {"O", 15.999405},
  {"S", 32.064787},
  {"Se", 78.963},
  {"Te", 127.603},
  {"Po", 209},
  {"F", 18.99840325},
  {"Cl", 35.452938},
  {"Br", 79.903528},
  {"I", 126.904473},
  {"At", 210},
  {"He", 4.0026022},
  {"Ne", 20.17976},
  {"Ar", 39.9481},
  {"Kr", 83.7982},
  {"Xe", 131.2936},
  {"Rn", 222},
  {"Sc", 44.9559126},
  {"Ti", 47.866749},
  {"V", 50.94151},
  {"Cr", 51.996133},
  {"Mn", 54.9380455},
  {"Fe", 55.845146},
  {"Co", 58.9331955},
  {"Ni", 58.69344},
  {"Cu", 63.5463},
  {"Zn", 63.5463},
  {"Y", 88.905852},
  {"Zr", 91.2242},
  {"Nb", 92.906382},
  {"Mo", 95.962},
  {"Tc", 98},
  {"Ru", 101.072},
  {"Rh", 102.905502},
  {"Pd", 106.421},
  {"Ag", 107.868151},
  {"Cd", 112.4118},
  {"La", 138.905477},
  {"Ce", 140.1161},
  {"Pr", 140.907652},
  {"Nd", 144.2423},
  {"Pm", 145},
  {"Sm", 150.362},
  {"Eu", 151.9641},
  {"Gd", 157.253},
  {"Tb", 158.925352},
  {"Dy", 162.5001},
  {"Ho", 164.930322},
  {"Er", 167.2593},
  {"Tm", 168.934212},
  {"Yb", 173.0545},
  {"Lu", 174.96681},
  {"Hf", 178.492},
  {"Ta", 180.947882},
  {"W", 183.841},
  {"Re", 186.2071},
  {"Os", 190.233},
  {"Ir", 192.2173},
  {"Pt", 195.0849},
  {"Au", 196.9665694},
  {"Hg", 200.592},
  {"Ac", 227},
  {"Th", 232.038062},
  {"Pa", 231.035882},
  {"U", 238.028913},
  {"Np", 237},
  {"Pu", 244},
  {"Am", 244},
  {"Cm", 247},
  {"Bk", 247},
  {"Cf", 251},
  {"Es", 252},
  {"Fm", 257},
  {"Md", 258},
  {"No", 259},
  {"Lr", 262},
  {"Rf", 267},
  {"Db", 268},
  {"Sg", 269},
  {"Bh", 270},
  {"Hs", 269},
  {"Cn", 285}
};

int RelAtomicMassTableSize()
{
  return sizeof(relamptable) / sizeof(relamptable[0]);
}

void getRelAtomicMass(int i, char **atom, double *relams)
{
  if(atom != NULL)
    (*atom) = relamptable[i].atom;

  if(relams != NULL)
    (*relams) = relamptable[i].atomval;
}

/*
 * Periodic Table organized with:
 * Atom name, covalent radii in angstrom from:
 *
 * Covalent radii revisited
 * Beatriz Cordero,  Verónica Gómez,  Ana E. Platero-Prats,  Marc Revés,
 * Jorge Echeverría,  Eduard Cremades,  Flavia Barragána and   Santiago Alvarez*a
 * Dalton Trans., 2008, 2832-2838
 */

PeriodicTable covradiustable[] = {
  {"H", 0.31},
  {"D", 2.01410178}, /*Deuterium isotope */
  {"Li", 1.28},
  {"Na", 1.66},
  {"K", 2.03},
  {"Rb", 2.20},
  {"Cs", 2.44},
  {"Fr", 2.60},
  {"Be", 0.96},
  {"Mg", 1.41},
  {"Ca", 1.76},
  {"Sr", 1.95},
  {"Ba", 2.15},
  {"Ra", 2.21},
  {"B", 0.84},
  {"Al", 1.21},
  {"Ga", 1.22},
  {"In", 1.42},
  {"Tl", 1.45},
  {"C.3", 0.76},
  {"C.2", 0.73},
  {"C.1", 0.69},
  {"Si", 1.11},
  {"Ge", 1.20},
  {"Sn", 1.39},
  {"Pb", 1.46},
  {"N", 0.71},
  {"P", 1.07},
  {"As", 1.19},
  {"Sb", 1.39},
  {"Bi", 1.48},
  {"O", 0.66},
  {"S", 1.05},
  {"Se", 1.20},
  {"Te", 1.38},
  {"Po", 1.40},
  {"F", 0.57},
  {"Cl", 1.02},
  {"Br", 1.20},
  {"I", 1.39},
  {"At", 1.50},
  {"He", 0.28},
  {"Ne", 0.58},
  {"Ar", 1.06},
  {"Kr", 1.16},
  {"Xe", 1.40},
  {"Rn", 1.50},
  {"Sc", 1.70},
  {"Ti", 1.60},
  {"V", 1.53},
  {"Cr", 1.39},
  {"Mn", 1.39}, /* Mn l.s.: 1.39; Mn h.s.: 1.61*/
  {"Fe", 1.32}, /*Fe l.s.: 1.32; Fe h.s: 1.52 */
  {"Co", 1.26}, /*Co l.s.: 1.26; Co h.s.: 1.50*/
  {"Ni", 1.24},
  {"Cu", 1.32},
  {"Zn", 1.22},
  {"Zr", 1.75},
  {"Nb", 1.64},
  {"Mo", 1.54},
  {"Tc", 1.47},
  {"Ru", 1.46},
  {"Rh", 1.42},
  {"Pd", 1.39},
  {"Ag", 1.45},
  {"Cd", 1.44},
  {"La", 2.07},
  {"Ce", 2.04},
  {"Pr", 2.03},
  {"Nd", 2.01},
  {"Pm", 1.99},
  {"Sm", 1.98},
  {"Eu", 1.98},
  {"Gd", 1.96},
  {"Tb", 1.94},
  {"Dy", 1.92},
  {"Ho", 1.92},
  {"Er", 1.89},
  {"Tm", 1.90},
  {"Y", 1.90},
  {"Yb", 1.87},
  {"Lu", 1.87},
  {"Hf", 1.75},
  {"Ta", 1.70},
  {"W", 1.62},
  {"Re", 1.51},
  {"Os", 1.44},
  {"Ir", 1.41},
  {"Pt", 1.36},
  {"Au", 1.36},
  {"Hg", 1.32},
  {"Ac", 2.15},
  {"Th", 2.06},
  {"Pa", 2.00},
  {"U", 1.96},
  {"Np", 1.90},
  {"Pu", 1.87},
  {"Am", 1.80},
  {"Cm", 1.69}
};

int CovRadiusTableSize()
{
  return sizeof(covradiustable) / sizeof(covradiustable[0]);
}

void getCovRadius(int i, char **atom, double *cradius)
{
  if(atom != NULL)
    (*atom) = covradiustable[i].atom;

  if(cradius != NULL)
    (*cradius) = covradiustable[i].atomval;
}

double getCovRadiusfromAtomName(char *_atom)
{
  char *atom = StrTrim(_atom);
  int i;
  const int sz = CovRadiusTableSize();

  /*if(strcmp(atom, "C") == 0){
    return relamptable[19].atomval;
  }
  else if(strcmp(atom, "H") == 0){
    return relamptable[0].atomval;
  }
  else if(strcmp(atom, "D") == 0){
    return relamptable[1].atomval;
  }
  else if(strcmp(atom, "O") == 0){
    return relamptable[29].atomval;
  }
  else if(strcmp(atom, "N") == 0){
    return relamptable[24].atomval;
  }
  else if(strcmp(atom, "S") == 0){
    return relamptable[30].atomval;
  }
  else if(strcmp(atom, "P") == 0){
    return relamptable[25].atomval;
  }
  else if(strcmp(atom, "F") == 0){
    return relamptable[34].atomval;
  }
  else if(strcmp(atom, "Cl") == 0){
    return relamptable[35].atomval;
  }
  else if(strcmp(atom, "Br") == 0){
    return relamptable[36].atomval;
  }
  else if(strcmp(atom, "I") == 0){
    return relamptable[37].atomval;
  }
  else if(strcmp(atom, "") == 0){
    return 0.;
  }
  else{*/
    int id = 0;
    for(i = 0; i < sz; i++){
      if(strcmp(atom, covradiustable[i].atom) == 0){
        id = i;
        break;
      }
    }
    return covradiustable[id].atomval;
  /*}*/
}

#define MAX_ATOMS 19
#define TYPEDIM 256

struct type{
  float eneg;
  char name[TYPEDIM];
};

/* Allen Electronegativity table
 *  expressed on the Pauling scale
 */
static struct type eneg_table[MAX_ATOMS] = {
  {2.300, "H"},
  {0.912, "Li"},
  {0.869, "Na"},
  {0.734, "K"},
  {1.576, "Be"},
  {1.293, "Mg"},
  {1.034, "Ca"},
  {2.051, "B"},
  {1.613, "Al"},
  {2.544, "C"},
  {2.916, "Si"},
  {3.066, "N"},
  {2.253, "P"},
  {3.610, "O"},
  {2.589, "S"},
  {4.193, "F"},
  {2.869, "Cl"},
  {2.685, "Br"},
  {2.359, "I"}
};

/*Get Atom Electronegativity*/
double getEnegativity(const char *type)
{
  int i;
  double eneg = -1.;
  for(i = 0; i < MAX_ATOMS; i++){
    if(strstr(type, eneg_table[i].name) != NULL){
      eneg = eneg_table[i].eneg;
      break;
    }
    else{
      continue;
    }
  }
  return eneg;
}

int getAtomNumberfromAtomName(char *_atom){
  char *atom = StrTrim(_atom);
  int i;
  const int sz = sizeof(relamptable) / sizeof(relamptable[0]);
  if(strcmp(atom, "C") == 0){
    return 12;
  }
  else if(strcmp(atom, "H") == 0){
    return 1;
  }
  else if(strcmp(atom, "D") == 0){
    return 1;
  }
  else if(strcmp(atom, "O") == 0){
    return 8;
  }
  else if(strcmp(atom, "N") == 0){
    return 7;
  }
  else if(strcmp(atom, "S") == 0){
    return 16;
  }
  else if(strcmp(atom, "P") == 0){
    return 15;
  }
  else if(strcmp(atom, "F") == 0){
    return 9;
  }
  else if(strcmp(atom, "Cl") == 0){
    return 17;
  }
  else if(strcmp(atom, "Br") == 0){
    return 35;
  }
  else if(strcmp(atom, "I") == 0){
    return 53;
  }
  else if(strcmp(atom, "") == 0){
    return 0.;
  }
  else{
    for(i = 0; i < sz; i++){
      if(strcmp(atom, relamptable[i].atom) == 0){
        return (int)floor(relamptable[i].atomval/2.f);
      }
    }
  }

  return -9999;
}


double getAtomWeightfromAtomName(char *_atom){
  char *atom = StrTrim(_atom);
  int i;
  const int sz = sizeof(relamptable) / sizeof(relamptable[0]);
  if(strcmp(atom, "C") == 0){
    return relamptable[19].atomval;
  }
  else if(strcmp(atom, "H") == 0){
    return relamptable[0].atomval;
  }
  else if(strcmp(atom, "D") == 0){
    return relamptable[1].atomval;
  }
  else if(strcmp(atom, "O") == 0){
    return relamptable[29].atomval;
  }
  else if(strcmp(atom, "N") == 0){
    return relamptable[24].atomval;
  }
  else if(strcmp(atom, "S") == 0){
    return relamptable[30].atomval;
  }
  else if(strcmp(atom, "P") == 0){
    return relamptable[25].atomval;
  }
  else if(strcmp(atom, "F") == 0){
    return relamptable[34].atomval;
  }
  else if(strcmp(atom, "Cl") == 0){
    return relamptable[35].atomval;
  }
  else if(strcmp(atom, "Br") == 0){
    return relamptable[36].atomval;
  }
  else if(strcmp(atom, "I") == 0){
    return relamptable[37].atomval;
  }
  else if(strcmp(atom, "") == 0){
    return 0.;
  }
  else{
    for(i = 0; i < sz; i++){
      if(strcmp(atom, relamptable[i].atom) == 0){
        return relamptable[i].atomval;
      }
    }
  }

  return -9999;
}

double getExactAtomWeightfromAtomName(char *_atom){
  char *atom = StrTrim(_atom);
  if(strcmp(atom, "C") == 0){
    return 12.;
  }
  else if(strcmp(atom, "H") == 0){
    return 1.00782503207;
  }
  else if(strcmp(atom, "D") == 0){
    return 2.0141017778;
  }
  else if(strcmp(atom, "Na") == 0){
    return 22.989769280929;
  }
  else if(strcmp(atom, "K") == 0){
    return 38.9637066820;
  }
  else if(strcmp(atom, "Ca") == 0){
    return 39.9625909822;
  }
  else if(strcmp(atom, "O") == 0){
    return  15.9949146195616;
  }
  else if(strcmp(atom, "N") == 0){
    return 14.00307400486;
  }
  else if(strcmp(atom, "S") == 0){
    return 31.9720710015;
  }
  else if(strcmp(atom, "P") == 0){
    return 30.9737616320;
  }
  else if(strcmp(atom, "F") == 0){
    return  18.998403227;
  }
  else if(strcmp(atom, "Cl") == 0){ /*Cl 35*/
    return  34.968852684;
  }
  else if(strcmp(atom, "Br") == 0){
    return 78.918337122;
  }
  else if(strcmp(atom, "I") == 0){
    return 126.9044734;
  }
  else if(strcmp(atom, "B") == 0){
    return 11.00930544;
  }
  else if(strcmp(atom, "Si") == 0){
    return 27.9769265325;
  }
  else if(strcmp(atom, "Ti") == 0){
    return 47.94794639;
  }
  else if(strcmp(atom, "Ni") == 0){
    return 57.93534297;
  }
  else if(strcmp(atom, "Se") == 0){
    return 79.916521321;
  }
  else if(strcmp(atom, "Sn") == 0){
    return 111.9048185;
  }
  else if(strcmp(atom, "Cu") == 0){
    return 62.92959756;
  }
  else if(strcmp(atom, "Fe") == 0){
    return 55.93493757;
  }
  else if(strcmp(atom, "Mg") == 0){
    return 23.98504170014;
  }
  else if(strcmp(atom, "Zn") == 0){
    return 63.92914227;
  }
  else if(strcmp(atom, "Mn") == 0){
    return 54.93804517;
  }
  else if(strcmp(atom, "Co") == 0){
    return 58.93319507;
  }
  else if(strcmp(atom, "Cr") == 0){
    return 51.94050758;
  }
  else if(strcmp(atom, "Li") == 0){
    return 7.016004558;
  }
  else if(strcmp(atom, "As") == 0){
    return 74.921596520;
  }
  else if(strcmp(atom, "Cd") == 0){
    return 113.903358529;
  }
  else if(strcmp(atom, "Hg") == 0){
    return 201.97064306;
  }
  else if(strcmp(atom, "Mo") == 0){
    return 97.905408221;
  }
  else if(strcmp(atom, "Ti") == 0){
    return 47.94794639;
  }
  else if(strcmp(atom, "V") == 0){
    return 50.943959511;
  }
  else if(strcmp(atom, "Al") == 0){
    return 26.9815386312;
  }
  else if(strcmp(atom, "Au") == 0){
    return 196.96656876;
  }
  else if(strcmp(atom, "Ba") == 0){
    return 137.90524725;
  }
  else if(strcmp(atom, "Pb") == 0){
    return 207.976652113;
  }
  else if(strcmp(atom, "W") == 0){
    return 183.95093129;
  }
  else if(strcmp(atom, "Be") == 0){
    return 9.01218224;
  }
  else if(strcmp(atom, "Ag") == 0){
    return 106.9050975;
  }
  else if(strcmp(atom, "Ge") == 0){
    return 73.921177818;
  }
  else if(strcmp(atom, "Rb") == 0){
    return 84.91178973812;
  }
  else if(strcmp(atom, "Rh") == 0){
    return 102.9055043;
  }
  else if(strcmp(atom, "Sb") == 0){
    return 120.903815724;
  }
  else if(strcmp(atom, "Ru") == 0){
    return 101.904349322;
  }
  else if(strcmp(atom, "Xe") == 0){
    return 131.90415509;
  }
  else if(strcmp(atom, "Ta") == 0){
    return 180.947995820;
  }
  else if(strcmp(atom, "Tl") == 0){
    return 204.974427514;
  }
  else if(strcmp(atom, "Ra") == 0){
    return 226.025409825;
  }
  else if(strcmp(atom, "Sr") == 0){
    return 87.905612257197;
  }
  else if(strcmp(atom, "Pt") == 0){
    return 194.96479119;
  }
  else if(strcmp(atom, "Te") == 0){
    return 129.906224421;
  }
  else if(strcmp(atom, "Bi") == 0){
    return 208.980398716;
  }
  else if(strcmp(atom, "Rn") == 0){
    return 208.980398716;
  }
  else if(strcmp(atom, "Gd") == 0){
    return 159.927054127;
  }
  else if(strcmp(atom, "") == 0){
    return 0.;
  }
  else{
    return -9999999;
  }
}

/*Van der Waals radii expressed in Angstrom  1 A = 1E-10 meters */
double getVanDerWaalsRadiifromAtomName(char *_atom){
  char *atom = StrTrim(_atom);
  if(strcmp(atom, "C") == 0){
    return 1.70;
  }
  else if(strcmp(atom, "H") == 0){
    return 1.2;
  }
  else if(strcmp(atom, "Li") == 0){
    return 1.82;
  }
  else if(strcmp(atom, "Na") == 0){
    return 2.77;
  }
  else if(strcmp(atom, "K") == 0){
    return 2.75;
  }
  else if(strcmp(atom, "Mg") == 0){
    return 1.73;
  }
  else if(strcmp(atom, "O") == 0){
    return 1.52;
  }
  else if(strcmp(atom, "N") == 0){
    return 1.55;
  }
  else if(strcmp(atom, "S") == 0){
    return 1.80;
  }
  else if(strcmp(atom, "P") == 0){
    return 1.80;
  }
  else if(strcmp(atom, "F") == 0){
    return 1.47;
  }
  else if(strcmp(atom, "Cl") == 0){
    return 1.75;
  }
  else if(strcmp(atom, "Br") == 0){
    return 1.85;
  }
  else if(strcmp(atom, "I") == 0){
    return 1.98;
  }
  else if(strcmp(atom, "Si") == 0){
    return 2.10;
  }
  else if(strcmp(atom, "Cu") == 0){
    return 1.40;
  }
  else if(strcmp(atom, "Ag") == 0){
    return 1.72;
  }
  else if(strcmp(atom, "Zn") == 0){
    return 1.39;
  }
  else if(strcmp(atom, "Ni") == 0){
    return 1.63;
  }
  else if(strcmp(atom, "Se") == 0){
    return 1.90;
  }
  else if(strcmp(atom, "") == 0){
    return 0.;
  }
  else{
    return -9999999;
  }
}


void getVanDerWaalsParams(char *_atom, double *D_ii, double *x_ii)
{
  /*
   * Parameters come from
   * doi:10.1.1.208.7677
   * UFF, a Full Periodic Table Force Field for
   * Molecular Mechanics and Molecular Dynamics Simulations
   * J. Am. Chem. SOC., Vol. 114, No. 25, 1992
   * N.B.: Combination rules:
   * x_ij = sqrt(x_ii * x_jj)
   * D_ij = sqrt(D_ii*D_jj)
   */
  char *atom = StrTrim(_atom);
  if(strcmp(atom, "C") == 0){
    (*D_ii) = 0.105;
    (*x_ii) = 3.851;
  }
  else if(strcmp(atom, "H") == 0){
    (*D_ii) = 0.044;
    (*x_ii) = 2.886;
  }
  else if(strcmp(atom, "Li") == 0){
    (*D_ii) = 0.025;
    (*x_ii) = 2.451;
  }
  else if(strcmp(atom, "Na") == 0){
    (*D_ii) = 0.030;
    (*x_ii) = 2.983;
  }
  else if(strcmp(atom, "K") == 0){
    (*D_ii) = 0.035;
    (*x_ii) = 3.812;
  }
  else if(strcmp(atom, "Mg") == 0){
    (*D_ii) = 0.111;
    (*x_ii) = 3.021;
  }
  else if(strcmp(atom, "O") == 0){
    (*D_ii) = 0.060;
    (*x_ii) = 3.500;
  }
  else if(strcmp(atom, "N") == 0){
    (*D_ii) = 0.069;
    (*x_ii) = 3.660;
  }
  else if(strcmp(atom, "S") == 0){
    (*D_ii) = 0.274;
    (*x_ii) = 4.035;
  }
  else if(strcmp(atom, "P") == 0){
    (*D_ii) = 0.305;
    (*x_ii) = 4.147;
  }
  else if(strcmp(atom, "F") == 0){
    (*D_ii) = 0.050;
    (*x_ii) = 3.364;
  }
  else if(strcmp(atom, "Cl") == 0){
    (*D_ii) = 0.227;
    (*x_ii) = 3.189;
  }
  else if(strcmp(atom, "Br") == 0){
    (*D_ii) = 0.251;
    (*x_ii) = 4.189;
  }
  else if(strcmp(atom, "I") == 0){
    (*D_ii) = 0.339;
    (*x_ii) = 4.50;
  }
  else if(strcmp(atom, "Si") == 0){
    (*D_ii) = 0.402;
    (*x_ii) = 4.295;
  }
  else if(strcmp(atom, "Au") == 0){
    (*D_ii) = 0.039;
    (*x_ii) = 3.293;
  }
  else if(strcmp(atom, "Al") == 0){
    (*D_ii) = 0.505;
    (*x_ii) = 4.499;
  }
  else if(strcmp(atom, "Ag") == 0){
    (*D_ii) = 0.036;
    (*x_ii) = 3.148;
  }
  else if(strcmp(atom, "Fe") == 0){
    (*D_ii) = 0.013;
    (*x_ii) = 2.912;
  }
  else if(strcmp(atom, "Cu") == 0){
    (*D_ii) = 0.005;
    (*x_ii) = 3.495;
  }
  else if(strcmp(atom, "Ag") == 0){
    (*D_ii) = 0.036;
    (*x_ii) = 3.148;
  }
  else if(strcmp(atom, "Zn") == 0){
    (*D_ii) = 0.124;
    (*x_ii) = 2.763;
  }
  else if(strcmp(atom, "Ni") == 0){
    (*D_ii) = 0.015;
    (*x_ii) = 2.834;
  }
  else if(strcmp(atom, "Se") == 0){
    (*D_ii) = 0.291;
    (*x_ii) = 4.205;
  }
  else if(strcmp(atom, "") == 0){
    (*D_ii) = 0.0;
    (*x_ii) = 0.0;
  }
  else{
    (*D_ii) = 0.0;
    (*x_ii) = 0.0;
  }
}


double getStaticScalarDipolePolarizability(char *_atom)
{
  /*
   * Parameters come from
   * doi: 10.1080/00268976.2018.1535143
   */
  char *atom = StrTrim(_atom);
  if(strcmp(atom, "C") == 0){
    return 11.3000;
  }
  else if(strcmp(atom, "H") == 0){
    return 4.50456;
  }
  else if(strcmp(atom, "Li") == 0){
    return 164.1125;
  }
  else if(strcmp(atom, "Be") == 0){
    return 37.74;
  }
  else if(strcmp(atom, "B") == 0){
    return 20.50;
  }
  else if(strcmp(atom, "Na") == 0){
    return 162.7;
  }
  else if(strcmp(atom, "K") == 0){
    return 289.7;
  }
  else if(strcmp(atom, "Mg") == 0){
    return 71.2;
  }
  else if(strcmp(atom, "O") == 0){
    return 5.3;
  }
  else if(strcmp(atom, "N") == 0){
    return 7.4;
  }
  else if(strcmp(atom, "S") == 0){
    return 19.4;
  }
  else if(strcmp(atom, "P") == 0){
    return 25;
  }
  else if(strcmp(atom, "F") == 0){
    return 3.74;
  }
  else if(strcmp(atom, "Cl") == 0){
    return 14.6;
  }
  else if(strcmp(atom, "Br") == 0){
    return 21;
  }
  else if(strcmp(atom, "I") == 0){
    return 32.9;
  }
  else if(strcmp(atom, "Si") == 0){
    return 37.3;
  }
  else if(strcmp(atom, "Ca") == 0){
    return 160.8;
  }
  else if(strcmp(atom, "Au") == 0){
    return 36;
  }
  else if(strcmp(atom, "Pb") == 0){
    return 47;
  }
  else if(strcmp(atom, "Al") == 0){
    return 57.8;
  }
  else if(strcmp(atom, "Ag") == 0){
    return 55;
  }
  else if(strcmp(atom, "Fe") == 0){
    return 62;
  }
  else if(strcmp(atom, "Sc") == 0){
    return 97;
  }
  else if(strcmp(atom, "Ti") == 0){
    return 100;
  }
  else if(strcmp(atom, "V") == 0){
    return 87;
  }
  else if(strcmp(atom, "Cr") == 0){
    return 83;
  }
  else if(strcmp(atom, "Mn") == 0){
    return 68;
  }
  else if(strcmp(atom, "Co") == 0){
    return 55;
  }
  else if(strcmp(atom, "Ni") == 0){
    return 49;
  }
  else if(strcmp(atom, "Cu") == 0){
    return 46.5;
  }
  else if(strcmp(atom, "Zn") == 0){
    return 38.67;
  }
  else if(strcmp(atom, "Ga") == 0){
    return 50;
  }
  else if(strcmp(atom, "Ge") == 0){
    return 40;
  }
  else if(strcmp(atom, "As") == 0){
    return 30;
  }
  else if(strcmp(atom, "Se") == 0){
    return 28.9;
  }
  else if(strcmp(atom, "Kr") == 0){
    return 16.78;
  }
  else if(strcmp(atom, "Rb") == 0){
    return 318.8;
  }
  else if(strcmp(atom, "Sr") == 0){
    return 197.2;
  }
  else if(strcmp(atom, "Y") == 0){
    return 162;
  }
  else if(strcmp(atom, "Zr") == 0){
    return 112;
  }
  else{
    return 0.0;
  }
}

double getFirstIonizationPotential(char *_atom)
{
  /*
   * Parameters come from
   * https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
   * and are expressed in eV
   */
  char *atom = StrTrim(_atom);
  if(strcmp(atom, "C") == 0){
    return 11.26030;
  }
  else if(strcmp(atom, "H") == 0){
    return 13.59844;
  }
  else if(strcmp(atom, "Li") == 0){
    return 5.39172;
  }
  else if(strcmp(atom, "Be") == 0){
    return 9.3227;
  }
  else if(strcmp(atom, "B") == 0){
    return 8.29803;
  }
  else if(strcmp(atom, "Na") == 0){
    return 5.13908;
  }
  else if(strcmp(atom, "K") == 0){
    return 4.34066;
  }
  else if(strcmp(atom, "Mg") == 0){
    return 7.64624;
  }
  else if(strcmp(atom, "O") == 0){
    return 13.61806;
  }
  else if(strcmp(atom, "N") == 0){
    return 14.53414;
  }
  else if(strcmp(atom, "S") == 0){
    return 10.36001;
  }
  else if(strcmp(atom, "P") == 0){
    return 10.48669;
  }
  else if(strcmp(atom, "F") == 0){
    return 17.42282;
  }
  else if(strcmp(atom, "Cl") == 0){
    return 12.96764;
  }
  else if(strcmp(atom, "Br") == 0){
    return 11.81381;
  }
  else if(strcmp(atom, "I") == 0){
    return 10.45126;
  }
  else if(strcmp(atom, "Si") == 0){
    return 8.15169;
  }
  else if(strcmp(atom, "Ca") == 0){
    return 6.11316;
  }
  else if(strcmp(atom, "Au") == 0){
    return 9.2255;
  }
  else if(strcmp(atom, "Pb") == 0){
    return 7.41666;
  }
  else if(strcmp(atom, "Al") == 0){
    return 5.98577;
  }
  else if(strcmp(atom, "Ag") == 0){
    return 7.5762;
  }
  else if(strcmp(atom, "Fe") == 0){
    return 7.9024;
  }
  else if(strcmp(atom, "Sc") == 0){
    return 6.5615;
  }
  else if(strcmp(atom, "Ti") == 0){
    return 6.8281;
  }
  else if(strcmp(atom, "V") == 0){
    return 6.7462;
  }
  else if(strcmp(atom, "Cr") == 0){
    return 6.7665;
  }
  else if(strcmp(atom, "Mn") == 0){
    return 7.43402;
  }
  else if(strcmp(atom, "Co") == 0){
    return 7.8810;
  }
  else if(strcmp(atom, "Ni") == 0){
    return 7.6398;
  }
  else if(strcmp(atom, "Cu") == 0){
    return 7.72638;
  }
  else if(strcmp(atom, "Zn") == 0){
    return 9.3942;
  }
  else if(strcmp(atom, "Ga") == 0){
    return 5.99930;
  }
  else if(strcmp(atom, "Ge") == 0){
    return 7.8994;
  }
  else if(strcmp(atom, "As") == 0){
    return 9.7886;
  }
  else if(strcmp(atom, "Se") == 0){
    return 9.75238;
  }
  else if(strcmp(atom, "Kr") == 0){
    return 13.99961;
  }
  else if(strcmp(atom, "Rb") == 0){
    return 4.17713;
  }
  else if(strcmp(atom, "Sr") == 0){
    return 5.6949;
  }
  else if(strcmp(atom, "Y") == 0){
    return 6.2171;
  }
  else if(strcmp(atom, "Zr") == 0){
    return 6.63390;
  }
  else{
    return 0.0;
  }
}


/* Atom Lonepairs
 * This function will return to you the atom lone pair
 * 1: one lonepairs
 * 2: two lonepairs
 * 3: three lonepairs
 */
int getAtomLonepairs(char *_atom, int connectivity, int bonds){
  //printf("%s %d %d\n", _atom, connectivity, bonds);
  char *atom = Trim(_atom);
  if(strcmp(atom, "O") == 0 && (connectivity == 2 || connectivity == 1) && (bonds == 2 || bonds == 0)){ // C=O or R-O-H or R-O-R or R-O-Ar or Ar-O-Ar
    return 2;
  }
  else if(strcmp(atom, "O") == 0 && connectivity == 2 && bonds == 3){ // Oxigen in aromatic ring with double bond
    return 1;
  }
  /* Connectivity 3: Amines
   * Connectivity 2: R-N=O
   * Connectivity 1: -CN (Nitrile)
   */
  else if(strcmp(atom, "N") == 0 && (connectivity == 3
                                    || connectivity == 2
                                    || connectivity == 1) && (bonds == 3 || bonds == 1 || bonds == 0)){ // N amines and temporary N in aromatic rings
    return 1;
  }
  else if(strcmp(atom, "N") == 0 && (connectivity == 4 || connectivity == 3) && bonds == 4){ // N amines
    return 0;
  }
  else if(strcmp(atom, "S") == 0 && (connectivity == 6 || connectivity == 4) && bonds == 6){ /* i.e.: SF6 D2SP3  and H2SO4 */
    return 0;
  }
  else if(strcmp(atom, "S") == 0 && (connectivity == 2 || connectivity == 1) && (bonds == 2 || bonds == 0)){ /* S=O, -S- or Ar-S-Ar*/
    return 2;
  }
  else if(strcmp(atom, "S") == 0 && connectivity == 3 && bonds == 4){ /* O-O-S=O */
    return 1;
  }
  else if(strcmp(atom, "P") == 0 && connectivity == 3 && (bonds == 3 || bonds == 1)){
    return 1;
  }
  else if(strcmp(atom, "P") == 0 && connectivity == 5 && bonds == 5){ /* i.e.: PCl5 D1SP3 */
    return 0;
  }
  else if((strcmp(atom, "F") == 0
    || strcmp(atom, "Cl") == 0
    || strcmp(atom, "Br") == 0
    || strcmp(atom, "I") == 0) && connectivity == 1 && bonds == 1){
    return 3;
  }
  else{
    return 0;
  }
}

int getAtomValenceElectrons(char *_atom){
  char *atom = Trim(_atom);
  if(strcmp(atom, "F") == 0 ||
           strcmp(atom, "Cl") == 0 ||
           strcmp(atom, "Br") == 0 ||
           strcmp(atom, "I") == 0 ||
           strcmp(atom, "As") == 0){
    return 7;
  }
  else if(strcmp(atom, "O") == 0 ||
     strcmp(atom, "S") == 0 ||
     strcmp(atom, "Se") == 0 ||
     strcmp(atom, "Te") == 0 ||
     strcmp(atom, "Po") == 0){
    return 6;
  }
  else if(strcmp(atom, "N") == 0 ||
          strcmp(atom, "P") == 0 ||
          strcmp(atom, "As") == 0 ||
          strcmp(atom, "Sb") == 0 ||
          strcmp(atom, "Bi") == 0){
    return 5;
  }
  else if(strcmp(atom, "C") == 0 ||
     strcmp(atom, "Si") == 0 ||
     strcmp(atom, "Ge") == 0 ||
     strcmp(atom, "Sn") == 0 ||
     strcmp(atom, "Pb") == 0){
    return 4;
  }
  else if(strcmp(atom, "B") == 0 ||
     strcmp(atom, "Al") == 0 ||
     strcmp(atom, "Ga") == 0 ||
     strcmp(atom, "In") == 0 ||
     strcmp(atom, "Tl") == 0){
    return 3;
  }
  else if(strcmp(atom, "Be") == 0 ||
     strcmp(atom, "Mg") == 0 ||
     strcmp(atom, "Ca") == 0 ||
     strcmp(atom, "Sr") == 0 ||
     strcmp(atom, "Ba") == 0 ||
     strcmp(atom, "Ra") == 0){
    return 2;
  }
  else if(strcmp(atom, "H") == 0 ||
     strcmp(atom, "Li") == 0 ||
     strcmp(atom, "Na") == 0 ||
     strcmp(atom, "K") == 0 ||
     strcmp(atom, "Rb") == 0 ||
     strcmp(atom, "Cs") == 0 ||
     strcmp(atom, "Fr") == 0){
    return 1;
  }
  else if(strcmp(atom, "Sc") == 0 ||
     strcmp(atom, "Y") == 0){
    return 3;
  }
  else if(strcmp(atom, "Ti") == 0 ||
          strcmp(atom, "Zr") == 0 ||
          strcmp(atom, "Hf") == 0){
    return 4;
  }
  else if(strcmp(atom, "V") == 0 ||
          strcmp(atom, "Nb") == 0 ||
          strcmp(atom, "Ta") == 0){
    return 5;
  }
  else if(strcmp(atom, "Cr") == 0 ||
          strcmp(atom, "Mo") == 0 ||
          strcmp(atom, "W") == 0){
    return 6;
  }
  else if(strcmp(atom, "Mn") == 0 ||
          strcmp(atom, "Tc") == 0 ||
          strcmp(atom, "Re") == 0){
    return 7;
  }
  else if(strcmp(atom, "Fe") == 0 ||
          strcmp(atom, "Ru") == 0 ||
          strcmp(atom, "Os") == 0){
    return 8;
  }
  else if(strcmp(atom, "Co") == 0 ||
          strcmp(atom, "Rh") == 0 ||
          strcmp(atom, "Ir") == 0){
    return 9;
  }
  else if(strcmp(atom, "Ni") == 0 ||
          strcmp(atom, "Pd") == 0 ||
          strcmp(atom, "Pt") == 0){
    return 10;
  }
  else if(strcmp(atom, "Cu") == 0 ||
          strcmp(atom, "Ag") == 0 ||
          strcmp(atom, "Au") == 0){
    return 11;
  }
  else if(strcmp(atom, "Zn") == 0 ||
          strcmp(atom, "Cd") == 0 ||
          strcmp(atom, "Hg") == 0){
    return 12;
  }
  else{
    return 0;
  }
}

void ReadAtomProperties(char *filename, AtomsProperty **lst)
{
  size_t i;
  char l[LINE_MAX];
  FILE *f;

  if((f = fopen (filename, "r")) == NULL){
    printf ("File %s not found!!!", filename);
    return;
  }

  i = 0;
  while(fgets(l, LINE_MAX, f) != NULL){
    if(strstr(l, "#") != NULL){
      continue;
    }
    else{
      i++;
    }
  }
  fclose(f);

  (*lst) = malloc(sizeof(AtomsProperty));
  (*lst)->atoms = malloc(sizeof(AProp)*i);
  (*lst)->natoms = i;

  if((f = fopen (filename, "r")) == NULL){
    printf ("File %s not found!!!", filename);
    return;
  }

  i = 0;
  while(fgets(l, LINE_MAX, f) != NULL){
    if(strstr(l, "#") != NULL){
      continue;
    }
    else{
      Trim(l);
      //printf("%s\n", l);
      sscanf(l, "%s %lf", (*lst)->atoms[i].type, &((*lst)->atoms[i].property));
      //printf("%s %e\n", (*lst)->atoms[i].type, (*lst)->atoms[i].property);
      i++;
    }
  }
  fclose(f);
}

void DeleteAtomProperties(AtomsProperty **lst)
{
  /*size_t i;
  for(i = 0; i < (*lst)->natoms; i++){
    free((*lst)->atoms[i].symbol);
  }*/
  free((*lst)->atoms);
  free((*lst));
}

double getGenericAProperty(char *_type, AtomsProperty *lst)
{
  size_t i;
  char *atom = Trim(_type);
  for(i = 0; i < lst->natoms; i++){
    if(strcmp(atom, lst->atoms[i].type) == 0){
      return lst->atoms[i].property;
    }
    else{
      continue;
    }
  }
  return 9999.f;
}

int isatom(char *_atom)
{
  char *atom = StrTrim(_atom);
  int i;
  const int sz = sizeof(relamptable) / sizeof(relamptable[0]);
  if(strcmp(atom, "C") == 0){
    return 1;
  }
  else if(strcmp(atom, "H") == 0){
    return 1;
  }
  else if(strcmp(atom, "D") == 0){
    return 1;
  }
  else if(strcmp(atom, "O") == 0){
    return 1;
  }
  else if(strcmp(atom, "N") == 0){
    return 1;
  }
  else if(strcmp(atom, "S") == 0){
    return 1;
  }
  else if(strcmp(atom, "P") == 0){
    return 1;
  }
  else if(strcmp(atom, "F") == 0){
    return 1;
  }
  else if(strcmp(atom, "Cl") == 0){
    return 1;
  }
  else if(strcmp(atom, "Br") == 0){
    return 1;
  }
  else if(strcmp(atom, "I") == 0){
    return 1;
  }
  else if(strcmp(atom, "") == 0){
    return 0;
  }
  else{
    for(i = 0; i < sz; i++){
      if(strcmp(atom, relamptable[i].atom) == 0){
        return 1;
      }
      else{
        continue;
      }
    }
  }
  return 0;
}

/*
 * Isotopic table
 */

Isotope isotopes_table[] = {
  {"C", 2, {12., 13.0033548378, 0., 0., 0., 0., 0., 0., 0.}, {98.938, 1.078, 0., 0., 0., 0., 0., 0., 0.}},
  {"H", 2, {1.00782503207, 2.0141017778, 0., 0., 0., 0., 0., 0., 0.}, {99.985, 0.015, 0., 0., 0., 0., 0., 0., 0.}}, // WARNIG: GET NEW VALUES FOR ISOTOPES
  {"N", 2, {14.00307400486, 15.000109, 0., 0., 0., 0., 0., 0., 0.}, {99.63, 0.37, 0., 0., 0., 0., 0., 0., 0.}}, // WARNIG: GET NEW VALUES FOR ISOTOPES
  {"O", 3, {15.9949146195616, 16.999131, 17.999159, 0., 0., 0., 0., 0., 0.}, {99.76, 0.038, 0.2, 0., 0., 0., 0., 0., 0.}}, // WARNIG: GET NEW VALUES FOR ISOTOPES
  {"S", 3, {31.9720710015, 32.971459, 33.967868, 0., 0., 0., 0., 0., 0.}, {95.02, 0.75, 4.21, 0., 0., 0., 0., 0., 0.}}, // WARNIG: GET NEW VALUES FOR ISOTOPES
  {"P", 1, {30.9737616320, 0., 0., 0., 0., 0., 0., 0., 0.}, {100., 0., 0., 0., 0., 0., 0., 0., 0.}} // WARNIG: GET NEW VALUES FOR ISOTOPES
  };

Isotope getIsotopefromAtomName(char *_atom){
  char *atom = StrTrim(_atom);
  int i;
  const int sz = sizeof(isotopes_table) / sizeof(isotopes_table[0]);
  if(strcmp(atom, "C") == 0){
    return isotopes_table[0];
  }
  else if(strcmp(atom, "H") == 0){
    return isotopes_table[1];
  }
  else if(strcmp(atom, "N") == 0){
    return isotopes_table[2];
  }
  else if(strcmp(atom, "O") == 0){
    return isotopes_table[3];
  }
  else if(strcmp(atom, "S") == 0){
    return isotopes_table[4];
  }
  else if(strcmp(atom, "P") == 0){
    return isotopes_table[5];
  }
  else{
    for(i = 0; i < sz; i++){
      if(strcmp(atom, isotopes_table[i].name) == 0){
        return isotopes_table[i];
      }
    }
  }
  return (Isotope){"None", 0, {0., 0., 0., 0., 0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0., 0., 0., 0., 0.}};
}
