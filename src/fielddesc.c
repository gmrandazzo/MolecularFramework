/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "fielddesc.h"
#include "geomdesc.h"
#include "molecule.h"
#include "misc.h"
#include "periodic_table.h"

#include <scientific.h>
#include <math.h>

/* Rotational invariant spherical field analysis
 * Algorithm:
 * 1) Calculate the center of mass
 * 2) Set n_prev = 1 n_next = n_prev + 1 and build two sphere one of radius c+n_prev*radius_step and one with radius c+n_next*radius_step
 * 3) Count the number of points cotained between the two sphere and add the number to the spectra.
 * 4) Continue the calculaiton until the molecule is inside the sphere.
 */
void SphericalFieldAnalysis(MOLECULE molecule, matrix *field, double radius_step, matrix **spectra)
{
  size_t i;
  double ccx, ccy, ccz, tot_w;
  double max_x, max_y, max_z, dmax;
  double n_next, n_prev;
  dvector *row;
  NewDVector(&row, 4);

  /* Calculate the center of mass */
  ccx = ccy = ccz = tot_w = 0.f;
  max_x = -9999.f;
  max_y = -9999.f;
  max_z = -9999.f;
  for(i = 0; i < molecule.n_atoms; i++){
    double aw = getAtomWeightfromAtomName(molecule.atoms[i].asymbl);
    tot_w += aw;
    ccx += molecule.atoms[i].coord.x*aw;
    ccy += molecule.atoms[i].coord.y*aw;
    ccz += molecule.atoms[i].coord.z*aw;
    if(fabs(molecule.atoms[i].coord.x) > max_x )
      max_x = fabs(molecule.atoms[i].coord.x);

    if(fabs(molecule.atoms[i].coord.y) > max_y)
      max_y = fabs(molecule.atoms[i].coord.y);

    if(fabs(molecule.atoms[i].coord.z) > max_z)
      max_z = fabs(molecule.atoms[i].coord.z);
  }

  ccx /= tot_w;
  ccy /= tot_w;
  ccz /= tot_w;

  dmax = GetDistance(ccx, ccy, ccz, max_x, max_y, max_z);

  /* Set n_prev = 1 n_next = n_prev + 1
   * and build two sphere one of radius c+n_prev*radius_step
   * and one with radius c+n_next*radius_step
   */
  n_prev = 1;
  n_next = n_prev+1;
  while(1){
    double r_prev = n_prev*radius_step;
    double r_next = n_next*radius_step;
    if(r_next > dmax)
      break;
    double values = 0.f;
    int npnt = 0;
    double min_val = 9999.f;
    double max_val = -9999.f;
    for(i = 0; i < field->row; i++){
      double pntx = field->data[i][0];
      double pnty = field->data[i][1];
      double pntz = field->data[i][2];

      double d = GetDistance(ccx, ccy, ccz, pntx, pnty, pntz);
      if(d >= r_prev && d <= r_next){
        values += field->data[i][3];

        if(field->data[i][3] < min_val)
          min_val = field->data[i][3];

        if(field->data[i][3] > max_val)
          max_val = field->data[i][3];

        npnt++;
      }
      else{
        continue;
      }
    }

    row->data[0] = n_next;

    if(values != 0.f){
      values /= (double)npnt;
      row->data[1] = values;
    }
    else{
      row->data[1] = 0.f;
    }

    if(min_val != 9999.f){
      row->data[2] = min_val;
    }
    else{
      row->data[2] = 0.f;
    }

    if(max_val != -9999.f){
      row->data[3] = max_val;
    }
    else{
      row->data[3] = 0.f;
    }

    MatrixAppendRow(spectra, row);
    n_prev++;
    n_next++;
  }
  DelDVector(&row);
}


/*
 * Calculate the average of the interaction points between min and max
 */
double AvgFieldValue(matrix *field, double min, double max)
{
  double _min = min, _max = max;
  size_t i, npnt;
  double val = 0.f;
  npnt = 0;
  for(i = 0; i < field->row; i++){
    if(field->data[i][3] > _min && field->data[i][3] < _max){
      val += field->data[i][3];
      npnt++;
    }
    else
      continue;
  }

  if(npnt > 0)
    return val/(double)npnt;
  else
    return 0.f;
}

/*
 * Calculate the frequency of the interaction points between min and max
 */
size_t FreqFieldValue(matrix *field, double min, double max)
{
  double _min = min, _max = max;
  size_t i, npnt;
  npnt = 0;
  for(i = 0; i < field->row; i++){
    if(field->data[i][3] > _min && field->data[i][3] < _max){
      npnt++;
    }
    else
      continue;
  }

  if(npnt > 0)
    return npnt;
  else
    return 0;
}


/*
 * Calculate the Baricenter of charges between min and max and the baricenter of the molecule
 */
double BaricenterFieldValue(MOLECULE molecule, matrix *field, double min, double max)
{
  size_t i, npnt;
  double ccx, ccy, ccz;
  double pntx, pnty, pntz;

  ccx = ccy = ccz = 0.f;
  for(i = 0; i < molecule.n_atoms; i++){
   ccx += molecule.atoms[i].coord.x;
   ccy += molecule.atoms[i].coord.y;
   ccz += molecule.atoms[i].coord.z;
  }

  ccx /= (double)molecule.n_atoms;
  ccy /= (double)molecule.n_atoms;
  ccz /= (double)molecule.n_atoms;

  pntx = pnty = pntz = 0.f;
  npnt = 0;
  for(i = 0; i < field->row; i++){
   if(field->data[i][3] > min && field->data[i][3] < max){
     pntx += field->data[i][0];
     pnty += field->data[i][1];
     pntz += field->data[i][2];
     npnt++;
   }
   else
     continue;
  }

  pntx /= (double)npnt;
  pnty /= (double)npnt;
  pntz /= (double)npnt;

  if(npnt > 0)
    return (double)GetDistance(ccx, ccy, ccz, pntx, pnty, pntz) * npnt;
  else
    return 0.f;
}

/*
 * Calculate the difference of the interaction points between two fields (A; B)
 */
double DiffFieldValue(matrix *field_a, matrix *field_b, double min, double max)
{
  size_t i, npnt_a, npnt_b;
  npnt_a = 0;
  npnt_b = 0;
  for(i = 0; i < field_a->row; i++){
    if(field_a->data[i][3] > min && field_a->data[i][3] < max){
      npnt_a++;
    }
    else
      continue;
  }

  for(i = 0; i < field_b->row; i++){
    if(field_b->data[i][3] > min && field_b->data[i][3] < max){
      npnt_b++;
    }
    else
      continue;
  }
  return (double) npnt_a - npnt_b;
}

/*
 * Calculate the ration of the interaction points between two fields
 */
double RatioFieldValue(matrix *field_a, matrix *field_b, double min, double max)
{
  size_t i, npnt_a, npnt_b;
  npnt_a = 0;
  npnt_b = 0;
  for(i = 0; i < field_a->row; i++){
   if(field_a->data[i][3] > min && field_a->data[i][3] < max){
     npnt_a++;
   }
   else
     continue;
  }

  for(i = 0; i < field_b->row; i++){
   if(field_b->data[i][3] > min && field_b->data[i][3] < max){
     npnt_b++;
   }
   else
     continue;
  }

  if(npnt_b > 0 && npnt_a > 0)
    return ((double)npnt_a)/((double)npnt_b);
  else
    return 0.f;
}

/*
 * Calculate the Dipole of energy between min and max
 */
double DipoleFieldValue(matrix *field, double min, double max)
{
 size_t i, npnt;
 double upntx, upnty, upntz;

 upntx = upnty = upntz = 0.f;
 npnt = 0;
 for(i = 0; i < field->row; i++){
   if(field->data[i][3] > min && field->data[i][3] < max){
     upntx += field->data[i][0]*field->data[i][3];
     upnty += field->data[i][1]*field->data[i][3];
     upntz += field->data[i][2]*field->data[i][3];
     npnt++;
   }
   else
     continue;
 }
 return sqrt(square(upntx)+square(upnty)+square(upntz));
}
