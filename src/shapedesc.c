/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <scientific.h>

#include "atomanalysis.h"
#include "geomdesc.h"
#include "shapedesc.h"
#include "molecule.h"
#include "misc.h"

void initShapePoints(SHAPEPNT **pnt)
{
  (*pnt) = malloc(sizeof(SHAPEPNT));
  initMatrix(&(*pnt)->surf);
  initMatrix(&(*pnt)->vol);
}

void deleteShapePoints(SHAPEPNT **pnt)
{
  DelMatrix(&(*pnt)->surf);
  DelMatrix(&(*pnt)->vol);
  free((*pnt));
}

void getShapePoints(MOLECULE molecule, SHAPEPNT **pnt, enum RADIUS_TYPE rtype)
{
  int i, j, k, l;
  double xmin, xmax, ymin, ymax, zmin, zmax, minradius, maxradius;
  double stepsize = 0.25, dx, dy, dz;

  dvector *cc;

  xmin = xmax = molecule.atoms[0].coord.x;
  ymin = ymax = molecule.atoms[0].coord.y;
  zmin = zmax = molecule.atoms[0].coord.z;

  initShapePoints(pnt);

  minradius = 9999.;
  maxradius = 0.;
  for(l = 0; l < molecule.n_atoms; l++){
    if(rtype == vanderwaals){
      molecule.atoms[l].radius = getVanDerWaalsRadiifromAtomName(molecule.atoms[l].asymbl);
    }
    else{
      molecule.atoms[l].radius = getCovRadiusfromAtomName(molecule.atoms[l].asymbl);
    }

    if(molecule.atoms[l].radius < minradius)
      minradius = molecule.atoms[l].radius;

    if(molecule.atoms[l].radius > maxradius)
      maxradius = molecule.atoms[l].radius;

    //printf("%f \n", molecule.atoms[l].radius);

    if(molecule.atoms[l].coord.x < xmin)
      xmin = molecule.atoms[l].coord.x;
    if(molecule.atoms[l].coord.x > xmax)
      xmax = molecule.atoms[l].coord.x;

    if(molecule.atoms[l].coord.y < ymin)
      ymin = molecule.atoms[l].coord.y;
    if(molecule.atoms[l].coord.y > ymax)
      ymax = molecule.atoms[l].coord.y;

    if(molecule.atoms[l].coord.z < zmin)
      zmin = molecule.atoms[l].coord.z;
    if(molecule.atoms[l].coord.z > zmax)
      zmax = molecule.atoms[l].coord.z;
  }

  double xspan = (xmax-xmin)/4.;
  double yspan = (ymax-ymin)/4.;
  double zspan = (zmax-zmin)/4.;

  xmin -= xspan;
  xmax += xspan;

  ymin -= yspan;
  ymax += yspan;

  zmin -= zspan;
  zmax += zspan;

  dx = xmin;
  dy = ymin;
  dz = zmin;
  int nxiter = ceil((xmax - xmin)/stepsize);
  int nyiter = ceil((ymax - ymin)/stepsize);
  int nziter =  ceil((zmax - zmin)/stepsize);
  /*printf("%f %f %f %f %f %f %d %d %d\n", xmin, xmax, ymin, ymax, zmin, zmax, nxiter, nyiter, nziter);*/
  NewDVector(&cc, 3);
  for(i = 0; i < nxiter; i++){
    printf("i: %d \n", i);
    cc->data[0] = dx + (stepsize/2.);
    dy = ymin;
    for(j = 0; j < nyiter; j++){
      cc->data[1] = dy + (stepsize/2.);
      dz = zmin;
      for(k = 0; k < nziter; k++){
        cc->data[2] = dz + (stepsize/2.);
        for(l = 0; l < molecule.n_atoms; l++){
          double dcatom = GetDistance(cc->data[0], cc->data[1], cc->data[2], molecule.atoms[l].coord.x, molecule.atoms[l].coord.y, molecule.atoms[l].coord.z);
          if(FLOAT_EQ(dcatom, molecule.atoms[l].radius, 1e-2)){
            MatrixAppendRow(&(*pnt)->surf, cc);
          }

          if(dcatom < molecule.atoms[l].radius || FLOAT_EQ(dcatom, molecule.atoms[l].radius, 1e-2)){
             MatrixAppendRow(&(*pnt)->vol, cc);
            break;
          }
          else{
            continue;
          }
        }
        dz += stepsize;
      }
      dy += stepsize;
    }
    dx += stepsize;
  }
  DelDVector(&cc);
}


void WriteMol2SurfPoints(SHAPEPNT *pnt, char *outfile)
{
  size_t i;
  FILE *fp;
  fp = fopen(outfile, "w");
  fprintf(fp, "@<TRIPOS>MOLECULE\n");
  fprintf(fp, "Volume Points\n");
  //fprintf(fp, "       0    0    0\n");
  fprintf(fp, "%5zu   0 0 0 0\n", pnt->surf->row);
  fprintf(fp, "SMALL\n");
  fprintf(fp, "NO_CHARGES\n\n");
  fprintf(fp, "@<TRIPOS>ATOM\n");
  for(i = 0; i < pnt->surf->row; i++){
    fprintf(fp, "%7zu ****    %10.2f%10.2f%10.2f Du\n", i+1, pnt->surf->data[i][0], pnt->surf->data[i][1], pnt->surf->data[i][2]);
  }
  fclose(fp);
}

void WriteMol2VolPoints(SHAPEPNT *pnt, char *outfile)
{
  size_t i;
  FILE *fp;
  fp = fopen(outfile, "w");
  fprintf(fp, "@<TRIPOS>MOLECULE\n");
  fprintf(fp, "Volume Points\n");
  //fprintf(fp, "       0    0    0\n");
  fprintf(fp, "%5zu   0 0 0 0\n", pnt->vol->row);
  fprintf(fp, "SMALL\n");
  fprintf(fp, "NO_CHARGES\n\n");
  fprintf(fp, "@<TRIPOS>ATOM\n");
  for(i = 0; i < pnt->vol->row; i++){
    fprintf(fp, "%7zu ****    %10.2f%10.2f%10.2f Du\n", i+1, pnt->vol->data[i][0], pnt->vol->data[i][1], pnt->vol->data[i][2]);
  }
  fclose(fp);
}


/*
 * The  Van der Waals molecular volume is calculated dividing the space in small
 * cube and summing these cubes if they are inside or not one atom of the molecule.
 */
void GetVDWMolVolSurf(MOLECULE molecule, double *volume, double *surface, double stepsize)
{
  int i, j, k, l;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double dx, dy, dz;
  double cubevol = stepsize*stepsize*stepsize; //volume cube
  int surfacepoints = 0;
  int allsurfacepoints = 0;
  (*surface) = 0.f;
  (*volume) = 0.f;

  xmin = xmax = molecule.atoms[0].coord.x;
  ymin = ymax = molecule.atoms[0].coord.y;
  zmin = zmax = molecule.atoms[0].coord.z;


  for(l = 0; l < molecule.n_atoms; l++){
    molecule.atoms[l].radius = getVanDerWaalsRadiifromAtomName(molecule.atoms[l].asymbl);
    (*surface) += 4*3.141592653589793*square(molecule.atoms[l].radius);

    if(molecule.atoms[l].coord.x < xmin)
      xmin = molecule.atoms[l].coord.x;
    if(molecule.atoms[l].coord.x > xmax)
      xmax = molecule.atoms[l].coord.x;

    if(molecule.atoms[l].coord.y < ymin)
      ymin = molecule.atoms[l].coord.y;
    if(molecule.atoms[l].coord.y > ymax)
      ymax = molecule.atoms[l].coord.y;

    if(molecule.atoms[l].coord.z < zmin)
      zmin = molecule.atoms[l].coord.z;
    if(molecule.atoms[l].coord.z > zmax)
      zmax = molecule.atoms[l].coord.z;
  }

  /*adding a factor of  20 % and a van der waals radius maximum of 3 */
  double scalefact = 0.2;
  double vdwmaxradius = 3;
  if(xmin < 0 && xmax > 0){
    xmin += xmin*scalefact - vdwmaxradius;
    xmax += xmax*scalefact + vdwmaxradius;
  }
  else if(xmin < 0 && xmax < 0){
    xmin += xmin*scalefact - vdwmaxradius;
    xmax += fabs(xmax*scalefact) + vdwmaxradius;
  }
  else{
    xmin += -1*(xmin*scalefact + vdwmaxradius);
    xmax += xmax*scalefact + vdwmaxradius;
  }

  if(ymin < 0 && ymax > 0){
    ymin += ymin*scalefact - vdwmaxradius;
    ymax += ymax*scalefact + vdwmaxradius;
  }
  else if(ymin < 0 && ymax < 0){
    ymin += ymin*scalefact - vdwmaxradius;
    ymax += fabs(ymax*scalefact) + vdwmaxradius;
  }
  else{
    ymin += -1*(ymin*scalefact + vdwmaxradius);
    ymax += ymax*scalefact + vdwmaxradius;
  }

  if(zmin < 0 && zmax > 0){
    zmin += zmin*scalefact - vdwmaxradius;
    zmax += zmax*scalefact + vdwmaxradius;
  }
  else if(zmin < 0 && zmax < 0){
    zmin += zmin*scalefact - vdwmaxradius;
    zmax += fabs(zmax*scalefact) + vdwmaxradius;
  }
  else{
    zmin += -1*(zmin*scalefact + vdwmaxradius);
    zmax += zmax*scalefact + vdwmaxradius;
  }

  dx = xmin;
  dy = ymin;
  dz = zmin;
  int nxiter = ceil((xmax - xmin)/stepsize);
  int nyiter = ceil((ymax - ymin)/stepsize);
  int nziter =  ceil((zmax - zmin)/stepsize);
  /*
  FILE *fp;
  fp = fopen("volumepoints.mol2", "w");
  fprintf(fp, "@<TRIPOS>MOLECULE\n");
  fprintf(fp, "Volume Points\n");
  fprintf(fp, "       0    0    0\n");
  fprintf(fp, "SMALL\n");
  fprintf(fp, "NO_CHARGES\n\n");
  fprintf(fp, "@<TRIPOS>ATOM\n");
  int q = 0;
  */
  /*printf("%f %f %f %f %f %f %d %d %d\n", xmin, xmax, ymin, ymax, zmin, zmax, nxiter, nyiter, nziter);*/
  for(i = 0; i < nxiter; i++){
    double ccx = dx + (stepsize/2.);
    dy = ymin;
    for(j = 0; j < nyiter; j++){
      double ccy = dy + (stepsize/2.);
      dz = zmin;
      for(k = 0; k < nziter; k++){
        double ccz = dz + (stepsize/2.);
        for(l = 0; l < molecule.n_atoms; l++){
          double dcatom = GetDistance(ccx, ccy, ccz, molecule.atoms[l].coord.x, molecule.atoms[l].coord.y, molecule.atoms[l].coord.z);
          if(dcatom < molecule.atoms[l].radius || FLOAT_EQ(dcatom, molecule.atoms[l].radius, 1e-2)){
            (*volume) += cubevol;
            if(FLOAT_EQ(dcatom, molecule.atoms[l].radius, 1e-2)){
              surfacepoints++;
             }
            /*fprintf(fp, "%7d ****    %10.2f%10.2f%10.2f Du\n", q, ccx, ccy, ccz);
            q++;*/
            break;
          }
          else{
            continue;
          }
        }

        for(l = 0; l < molecule.n_atoms; l++){
          double dcatom = GetDistance(ccx, ccy, ccz, molecule.atoms[l].coord.x, molecule.atoms[l].coord.y, molecule.atoms[l].coord.z);
          if(FLOAT_EQ(dcatom, molecule.atoms[l].radius, 1e-2)){
              allsurfacepoints++;
          }
          else{
            continue;
          }
        }

        dz += stepsize;
      }
      dy += stepsize;
    }
    dx += stepsize;
  }
  /*fclose(fp);*/
  (*volume) = (*volume)*0.997 + 0.1135;
  if(surfacepoints != allsurfacepoints)
    (*surface) *= (double)surfacepoints/(double)allsurfacepoints;
}
