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
 * https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii
 * and organized as:
 * Atom name, Relative atomic mass
 */
PeriodicTable relamptable[] = {
  {"H", 1.007940754056},
  {"He", 4.002601932121},
  {"Li", 6.940036602916},
  {"Be", 9.012183065000},
  {"B", 10.811028046410},
  {"C", 12.010735896735},
  {"N", 14.006703211446},
  {"O", 15.999404924318},
  {"F", 18.998403162730},
  {"Ne", 20.180046380520},
  {"Na", 22.989769282000},
  {"Mg", 24.305051619837},
  {"Al", 26.981538530000},
  {"Si", 28.085498705706},
  {"P", 30.973761998420},
  {"S", 32.064787406127},
  {"Cl", 35.452937582608},
  {"Ar", 39.947798563582},
  {"K", 39.098300910086},
  {"Ca", 40.078022511018},
  {"Sc", 44.955908280000},
  {"Ti", 47.866744962722},
  {"V", 50.941465037425},
  {"Cr", 51.996131755434},
  {"Mn", 54.938043910000},
  {"Fe", 55.845144433866},
  {"Co", 58.933194290000},
  {"Ni", 58.693347109948},
  {"Cu", 63.546039945830},
  {"Zn", 65.377782529525},
  {"Ga", 69.723066072594},
  {"Ge", 72.627550164687},
  {"As", 74.921594570000},
  {"Se", 78.959388557014},
  {"Br", 79.903527780510},
  {"Kr", 83.797999995326},
  {"Rb", 85.467663595620},
  {"Sr", 87.616644469620},
  {"Y", 88.905840300000},
  {"Zr", 91.223641597060},
  {"Nb", 92.906373000000},
  {"Mo", 95.959788541188},
  {"Ru", 101.064940139160},
  {"Rh", 102.905498000000},
  {"Pd", 106.415327507340},
  {"Ag", 107.868149634557},
  {"Cd", 112.411557818268},
  {"In", 114.818086629446},
  {"Sn", 118.710112593011},
  {"Sb", 121.759783673480},
  {"Te", 127.603126484660},
  {"I", 126.904471900000},
  {"Xe", 131.292761447791},
  {"Cs", 132.905451961000},
  {"Ba", 137.326891628632},
  {"La", 138.905468873713},
  {"Ce", 140.115730737855},
  {"Pr", 140.907657600000},
  {"Nd", 144.241596031827},
  {"Sm", 150.366355711930},
  {"Eu", 151.964378126380},
  {"Gd", 157.252130646880},
  {"Tb", 158.925354700000},
  {"Dy", 162.499472819424},
  {"Ho", 164.930328800000},
  {"Er", 167.259082649669},
  {"Tm", 168.934217900000},
  {"Yb", 173.054150166317},
  {"Lu", 174.966814957855},
  {"Hf", 178.484978723400},
  {"Ta", 180.947875636227},
  {"W", 183.841777550513},
  {"Re", 186.206704545600},
  {"Os", 190.224859628240},
  {"Ir", 192.216051652100},
  {"Pt", 195.084456864931},
  {"Au", 196.966568790000},
  {"Hg", 200.599167034556},
  {"Tl", 204.383412839360},
  {"Pb", 207.216908063000},
  {"Bi", 208.980399100000},
  {"Th", 232.038055800000},
  {"Pa", 231.035884200000},
  {"Pm", 145},
  {"Tc", 98},
  {"Po", 209},
  {"At", 210},
  {"Rn", 222},
  {"Fr", 223},
  {"Ra", 226},
  {"Ac", 227},
  {"Np", 237},
  {"Pu", 244}
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

    int id = 0;
    for(i = 0; i < sz; i++){
      if(strcmp(atom, covradiustable[i].atom) == 0){
        id = i;
        break;
      }
    }
    return covradiustable[id].atomval;
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
    return relamptable[5].atomval;
  }
  else if(strcmp(atom, "H") == 0){
    return relamptable[0].atomval;
  }
  else if(strcmp(atom, "D") == 0){
    return 2.01410177812;
  }
  else if(strcmp(atom, "O") == 0){
    return relamptable[7].atomval;
  }
  else if(strcmp(atom, "N") == 0){
    return relamptable[6].atomval;
  }
  else if(strcmp(atom, "S") == 0){
    return relamptable[15].atomval;
  }
  else if(strcmp(atom, "P") == 0){
    return relamptable[14].atomval;
  }
  else if(strcmp(atom, "F") == 0){
    return relamptable[8].atomval;
  }
  else if(strcmp(atom, "Cl") == 0){
    return relamptable[16].atomval;
  }
  else if(strcmp(atom, "Br") == 0){
    return relamptable[34].atomval;
  }
  else if(strcmp(atom, "I") == 0){
    return relamptable[51].atomval;
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
    return 12.000000000000;
  }
  else if(strcmp(atom, "H") == 0){
    return 1.007825032230;
  }
  else if(strcmp(atom, "D") == 0){
    return 2.01410177812;
  }
  else if(strcmp(atom, "Na") == 0){
    return 22.989769282000;
  }
  else if(strcmp(atom, "K") == 0){
    return 38.963706486400;
  }
  else if(strcmp(atom, "Ca") == 0){
    return 39.962590863000;
  }
  else if(strcmp(atom, "O") == 0){
    return  15.994914619570;
  }
  else if(strcmp(atom, "N") == 0){
    return 14.003074004430;
  }
  else if(strcmp(atom, "S") == 0){
    return 31.972071174400;
  }
  else if(strcmp(atom, "P") == 0){
    return 30.973761998420;
  }
  else if(strcmp(atom, "F") == 0){
    return 18.998403162730;
  }
  else if(strcmp(atom, "Cl") == 0){ /*Cl 35*/
    return 34.968852682000;
  }
  else if(strcmp(atom, "Br") == 0){
    return 78.918337600000;
  }
  else if(strcmp(atom, "I") == 0){
    return 126.904471900000;
  }
  else if(strcmp(atom, "B") == 0){
    return 11.009305360000;
  }
  else if(strcmp(atom, "Si") == 0){
    return 27.976926534650;
  }
  else if(strcmp(atom, "Ti") == 0){
    return 47.947941980000;
  }
  else if(strcmp(atom, "Ni") == 0){
    return 57.935342410000;
  }
  else if(strcmp(atom, "Se") == 0){
    return 79.916521800000;
  }
  else if(strcmp(atom, "Sn") == 0){
    return 119.902201630000;
  }
  else if(strcmp(atom, "Cu") == 0){
    return 62.929597720000;
  }
  else if(strcmp(atom, "Fe") == 0){
    return 55.934936330000;
  }
  else if(strcmp(atom, "Mg") == 0){
    return 23.985041697000;
  }
  else if(strcmp(atom, "Zn") == 0){
    return 63.929142010000;
  }
  else if(strcmp(atom, "Mn") == 0){
    return 54.938043910000;
  }
  else if(strcmp(atom, "Co") == 0){
    return 58.933194290000;
  }
  else if(strcmp(atom, "Cr") == 0){
    return 51.940506230000;
  }
  else if(strcmp(atom, "Li") == 0){
    return 7.016003436600;
  }
  else if(strcmp(atom, "As") == 0){
    return 74.921594570000;
  }
  else if(strcmp(atom, "Cd") == 0){
    return 113.903365090000;
  }
  else if(strcmp(atom, "Hg") == 0){
    return 201.970643400000;
  }
  else if(strcmp(atom, "Mo") == 0){
    return 97.905404820000;
  }
  else if(strcmp(atom, "V") == 0){
    return 50.943957040000;
  }
  else if(strcmp(atom, "Al") == 0){
    return 26.981538530000;
  }
  else if(strcmp(atom, "Au") == 0){
    return 196.966568790000;
  }
  else if(strcmp(atom, "Ba") == 0){
    return 137.905247000000;
  }
  else if(strcmp(atom, "Pb") == 0){
    return 207.976652500000;
  }
  else if(strcmp(atom, "W") == 0){
    return 183.950930920000;
  }
  else if(strcmp(atom, "Be") == 0){
    return 9.012183065000;
  }
  else if(strcmp(atom, "Ag") == 0){
    return 106.905091600000;
  }
  else if(strcmp(atom, "Ge") == 0){
    return 73.921177761000;
  }
  else if(strcmp(atom, "Rb") == 0){
    return 84.911789737900;
  }
  else if(strcmp(atom, "Rh") == 0){
    return 102.905498000000;
  }
  else if(strcmp(atom, "Sb") == 0){
    return 120.903812000000;
  }
  else if(strcmp(atom, "Ru") == 0){
    return 101.904344100000;
  }
  else if(strcmp(atom, "Xe") == 0){
    return 131.904155085600;
  }
  else if(strcmp(atom, "Ta") == 0){
    return 180.947995800000;
  }
  else if(strcmp(atom, "Tl") == 0){
    return 204.974427800000;
  }
  else if(strcmp(atom, "Ra") == 0){
    return 226.025409825;
  }
  else if(strcmp(atom, "Sr") == 0){
    return 87.905612500000;
  }
  else if(strcmp(atom, "Pt") == 0){
    return 194.964791700000;
  }
  else if(strcmp(atom, "Te") == 0){
    return 129.906222748000;
  }
  else if(strcmp(atom, "Bi") == 0){
    return 208.980399100000;
  }
  else if(strcmp(atom, "Rn") == 0){
    return 208.980398716;
  }
  else if(strcmp(atom, "Gd") == 0){
    return 157.924112300000;
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
