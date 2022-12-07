/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "mol3Dfields.h"
#include "atomanalysis.h"
#include "geomdesc.h"
#include "miscdesc.h"
#include "molecule.h"
#include "periodic_table.h"
#include "atomdesc.h"
#include "forcefield.h"

#include <math.h>
#include <scientific.h>

#include <unistd.h>

#define C_e_ -1.6021765314E-19

void NewVoxel(VOXEL **v, size_t nx, size_t ny, size_t nz)
{
  size_t i, j, k;
  (*v) = malloc(sizeof(VOXEL));
  (*v)->nx = nx;
  (*v)->ny = ny;
  (*v)->nz = nz;
  (*v)->xmin = (*v)->xmax = (*v)->ymin = (*v)->ymax = (*v)->zmin = (*v)->zmax = 0.f;
  (*v)->pnt = malloc(sizeof(double**)*nx);
  for(i = 0; i < nx; i++){
    (*v)->pnt[i] = malloc(sizeof(double*)*ny);
    for(j = 0; j < ny; j++){
      (*v)->pnt[i][j] = malloc(sizeof(double)*nz);
      for(k = 0; k < nz; k++){
        (*v)->pnt[i][j][k] = 0.f;
      }
    }
  }
}

void DelVoxel(VOXEL **v)
{
  size_t i, j;
  for(i = 0; i < (*v)->nx; i++){
    for(j = 0; j < (*v)->ny; j++){
      free((*v)->pnt[i][j]);
    }
    free((*v)->pnt[i]);
  }
  free((*v)->pnt);
  free((*v));
}

void Voxel2Matrix(VOXEL *v, matrix *m)
{
  size_t i, j, k, c;
  ResizeMatrix(m, v->nx*v->ny*v->nz, 4);
  double x, dx, y, dy, z, dz;
  x = v->xmin;
  dx = (v->xmax-v->xmin)/(double)v->nx;
  dy = (v->ymax-v->ymin)/(double)v->ny;
  dz = (v->zmax-v->zmin)/(double)v->nz;
  c = 0;
  for(i = 0; i < v->nx; i++){
    y = v->ymin;
    for(j = 0; j < v->ny; j++){
      z = v->zmin;
      for(k = 0; k < v->nz; k++){
        m->data[c][0] = x;
        m->data[c][1] = y;
        m->data[c][2] = z;
        m->data[c][3] = v->pnt[i][j][k];
        c++;
        z += dz;
      }
      y += dy;
    }
    x += dx;
  }
}

double irand(int min, int max){
  return ((double)rand() / ((double)RAND_MAX + 1.0)) * (max - min) + min;
}

void HBSwitchFunction(double r, double angle, double *s_r, double *s_angle)
{
  /* r < 4.0 angstrom HB ON
   * r > 4.5 angstrom HB OFF
   * 4.0 A < r < 4.5 A --> Formula according CHARMM forcefield
    */
  if(r < 4.0){ /* ON */
    (*s_r) = 1.f;
  }
  else if(r > 4.5){ /* OFF */
    (*s_r) = 0.f;
  }
  else{
    /* formula according CHARMM forcefield */
    (*s_r) = (square(4.5-r)*(4.5+(2*r)-(3*4.0))) / 0.125; /*0.125 = (4.5-4.0)^3*/
  }

   /* r < 60 degree HB ON
    * r > 75 adegree HB OFF
    * 60 A < r < 75  A --> Formula according CHARMM forcefield
    * Limit according to this:
    * Distribution of hydrogen bond angles in molecular crystals
    * MASAMI HASEGAWA & HARUHIKO NODA
    * Letters to Nature
    * Nature 254, 212 (20 March 1975)
    * doi:10.1038/254212a0
    */

  if(angle < 60 ){ /* ON */
    (*s_angle) = 1.f;
  }
  else if(angle > 75){ /* OFF */
    (*s_angle) = 0.f;
  }
  else{
    /* formula */
    (*s_angle) = (square(75-angle)*(75+(2*angle)-(3*60))) / 328509.0; /*328509.0 = (75-60)^3*/
  }
}

/*
 * Calculate the Hydrogen Bond interaction energy in Joule
 * The probe is the acceptor and the molecule is the donor
 */
double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, POINT acceptor_coord, double acc_charge)
{
  /*printf("HBEnergyCalc\n");
  printf("charge: %f  dipole: %f\n", acc_charge, u);*/
  /*center of the dipole */
  double ucx = (h_atom.x + hetatm.x) / 2.f;
  double ucy = (h_atom.y + hetatm.y) / 2.f;
  double ucz = (h_atom.z + hetatm.z) / 2.f;

  /*Angle calculation
  POINT p1; p1.x = hetatm.x; p1.y = hetatm.y; p1.z = hetatm.z;
  POINT p2; p2.x = acceptor_coord.x; p2.y = acceptor_coord.y; p2.z = acceptor_coord.z;
  POINT c; c.x = h_atom.x; c.y = h_atom.y; c.z = h_atom.z;
  double angle = Rad2Grad(GetAngle(p1, p2, c));*/

  double angle = Rad2Grad(GetAngle(hetatm, acceptor_coord, h_atom));
  double dgpatom = GetDistance(acceptor_coord.x, acceptor_coord.y, acceptor_coord.z, ucx, ucy, ucz);

  /*switch functions */
  double s_dst = 0.f;
  double s_angle = 0.f;
  HBSwitchFunction(dgpatom, angle, &s_dst, &s_angle);
  double hb_energy = (au2C(acc_charge) * D2Cm(u) *cos(angle)/(4*_pi_*epsilon0_C_J_m*square(dgpatom*1E-10)));
  /*if(s_dst > 0 && s_angle > 0)
    printf("HB angle: %f; Distance: %f; S(theta): %f; S(dst): %f E(HB): %f\n", angle, dgpatom, s_angle, s_dst, hb_energy);*/
  return hb_energy*s_dst*s_angle;
}

/* Build coordinate for Hydrogen when the probe is an HBD */
double HBAEnergyCalc(char *patomtype, POINT probecoord, double u, POINT mol_coord, double mol_charge)
{
  /*printf("HBAEnergyCalc\n");*/
  size_t i;
  POINT hatom;
  double hb_pot = 0.f;
  for(i = 0; i < 3; i++){
    hatom.x = probecoord.x;
    hatom.y = probecoord.y;
    hatom.z = probecoord.z;
    /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
    if(strcmp(patomtype, "N") == 0){
      /*N-H distance 1.01 Angstrom */
      if(i == 0){
        hatom.x += 1.01;
      }
      else if(i == 1){
        hatom.y += 1.01;
      }
      else{
        hatom.z += 1.01;
      }
    }
    else if(strcmp(patomtype, "O") == 0){
      /*O-H distance 0.96 Angstrom */
      if(i == 0){
        hatom.x += 0.96;
      }
      else if(i == 1){
        hatom.y += 0.96;
      }
      else{
        hatom.z += 0.96;
      }
    }
    else if(strcmp(patomtype, "S") == 0){
      /*S-H distance 1.33 Angstrom */
      if(i == 0){
        hatom.x += 1.33;
      }
      else if(i == 1){
        hatom.y += 1.33;
      }
      else{
        hatom.z += 1.33;
      }
    }
    else if(strcmp(patomtype, "F") == 0){
      /*C-F distance 1.35 Angstrom */
      if(i == 0){
        hatom.x -= 1.35;
      }
      else if(i == 1){
        hatom.y -= 1.35;
      }
      else{
        hatom.z -= 1.35;
      }
    }
    else if(strcmp(patomtype, "Cl") == 0){
      /*C-Cl distance 1.70 Angstrom*/
      if(i == 0){
        hatom.x -= 1.70;
      }
      else if(i == 1){
        hatom.y -= 1.70;
      }
      else{
        hatom.z -= 1.70;
      }
    }
    else if(strcmp(patomtype, "Si") == 0){
      /* Si-H distance 1.48 Angstrom */
      if(i == 0){
        hatom.x += 1.48;
      }
      else if(i == 1){
        hatom.y += 1.48;
      }
      else{
        hatom.z += 1.48;
      }
    }
    else{
      /* Add simply 1.5 angstrom; */
      if(i == 0){
        hatom.x += 1.50;
      }
      else if(i == 1){
        hatom.y += 1.50;
      }
      else{
        hatom.z += 1.50;
      }
    }
    /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, POINT acceptor_coord, double acc_charge)*/
    hb_pot += HBEnergyCalc(probecoord, hatom, u, mol_coord, mol_charge);
  }
  return hb_pot /= 3.f;
}

int HBDCheck(char *atomtype)
{
  if(strcmp(atomtype, "N") == 0 ||
     strcmp(atomtype, "O") == 0 ||
     strcmp(atomtype, "S") == 0 ||
     strcmp(atomtype, "F") == 0)
    return 0;
  else
    return 1;
}

int HBACheck(char *atomtype){
  if(strcmp(atomtype, "O") == 0 ||
         strcmp(atomtype, "N") == 0 ||
         strcmp(atomtype, "S") == 0 ||
         strcmp(atomtype, "F") == 0 ||
         strcmp(atomtype, "Cl") == 0 ||
         strcmp(atomtype, "P") == 0 ||
         strcmp(atomtype, "Si") == 0)
    return 0;
  else
    return 1;
}

/* Perform an uniform sampling in radiants of the angles:
 * -theta: z angle respect to the xy plane between (0, pi)
 * -phi: angle respect x and y between (0, 2pi)
*/
void SphereCordaAngleGenerator(size_t nphi, size_t ntheta, matrix *angles)
{
  double phi, theta, d_phi, d_theta;
  d_phi = _pi_/(double)nphi;
  d_theta = (2*_pi_)/(double)ntheta;
  printf("%f %f\n", d_phi, d_theta);
  dvector *row;
  NewDVector(&row, 2);
  theta = 0.f;
  while(theta <= 2*_pi_){
    phi = 0.f;
    while(phi <= _pi_){
      row->data[0] = theta;
      row->data[1] = phi;
      phi += d_phi;
      MatrixAppendRow(angles, row);
    }
    theta += d_theta;
  }
  DelDVector(&row);
}

void GenEquidistributedPointsInSphere(int npnt, double r, matrix *pnt_in_sphere)
{
  /* Algorith:
     How to generate equidistributed points on the surface of a sphere
     Markus Deserno

     Ncount = 0
     alpha = 4*pi*r^2/N
     d = sqrt(alpha)
     Mtheta = round(pi/N)
     d_theta = pi/Mtheta
     d_phi = alpha/d_theta
     for each m in 0... Mtheta-1 do:
       theta = pi*(m+0.5)/Mtheta
       Mphi = round(2*pi*n/d_phi)
       for each n in Mphi-1 do:
         phi = 2*pi*n/Mphi
         x =  r*sin(theta)*cos(phi)
         y = r*sin(theta)*sin(phi)
         z = r*cos(theta)

  */
  dvector *row;
  NewDVector(&row, 3);
  int i, j;
  int N = npnt;
  double alpha = 4*_pi_*r*r/N;
  double d = sqrt(alpha);
  int Mtheta = round(_pi_/d);
  double d_theta = _pi_/Mtheta;
  double d_phi = alpha/d_theta;
  for(i = 0; i < Mtheta; i++){
    double theta = _pi_*(i+0.5)/(double)Mtheta;
    double Mphi = round(2.f*_pi_*sin(theta)/d_phi);
    for(j = 0; j < Mphi; j++){
      double phi = 2.f*_pi_*j/(double)Mphi;
      row->data[0] = r*sin(theta)*cos(phi);
      row->data[1] = r*sin(theta)*sin(phi);
      row->data[2] = r*cos(theta);
      MatrixAppendRow(pnt_in_sphere, row);
    }
  }
  DelDVector(&row);
}

/*
 * Calculate the interaction fields according the VDW radius sphere
 */
void FieldCalculator(ForceField ff,
                     MOLECULE *molecule,
                     int formal_charge,
                     size_t npnt,
                     enum RADIUS_TYPE rtype,
                     int probe_id,
                     double pdistance,
                     matrix *field){
  size_t i, j, k;
  double pradius = 0.f;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  //matrix *angles;
  matrix *pnts;

  NewDVector(&row, 4);

  //initMatrix(&angles);
  /* angles theta and phi utilised to perform an uniform surface point sampling */
  //SphereCordaAngleGenerator(npnt, npnt, &angles);

  /*Before computing interaction energies perform an atom analysis */
  AssignParams(molecule, ff, rtype, formal_charge);

  /* Set for each atom the choosed radius
  for(i = 0; i < molecule->n_atoms; i++){
   if(rtype == vanderwaals){
     molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
   }
   else{
     molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }
  */

  char *patomtype = strdup(ff.params[probe_id].atom_name);
  if(rtype == vanderwaals){
    pradius = getVanDerWaalsRadiifromAtomName(patomtype);
  }
  else{
    pradius = getCovRadiusfromAtomName(patomtype);
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);

    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */
      double x_pnt = pnts->data[i][0];
      double y_pnt = pnts->data[i][1];
      double z_pnt = pnts->data[i][2];

      row->data[0] = x_pnt + molecule->atoms[j].coord.x;
      row->data[1] = y_pnt + molecule->atoms[j].coord.y;
      row->data[2] = z_pnt + molecule->atoms[j].coord.z;

      //printf("%f %f %f\n", row->data[0],row->data[1],row->data[2]);
      // check if is a surface point or not....
      // get == 0: surface points
      int get = 0;
      for(k = 0; k < molecule->n_atoms; k++){
        if(k != j){
          double dgpatom = GetDistance(row->data[0], row->data[1], row->data[2],
                                       molecule->atoms[k].coord.x, molecule->atoms[k].coord.y, molecule->atoms[k].coord.z);
          if (dgpatom < molecule->atoms[k].radius){ // this point is penetrating an other atom so is not in the surface
            get = -1;
            break;
          }
          else{
            continue;
          }
        }
        else{
          continue;
        }
      }

      if(get == 0){
        /* This is a surface point.
         * We can calculate the potential of this point function of other atoms cotributions!
         * N.B.: Charges must be in atomic unit (au)!
         * 1 hartree (Eh) = 27.2113834 eV
         * 1 Hartree = 627.503 kcal/mol
         * 1 eV = 23.0605 kcal/mol
         * 1 eV = 1 J/C
         * 1 Debye = 3.336·10^(-30) C·m
         */

        /*electrostatic energy due to the charge-charge interaction*/
        double el_pot = 0.f;

        /* Lenard-Jones potential calcualted according to the classical equation */
        double lj_pot = 0.f;

        /*
         * When calculate the energy of interaction we will take into account all the
         * other atoms charges
         */
        for(k = 0; k < molecule->n_atoms; k++){
          /* we need to add the atom probe radius to simulate the contact sphere system */
          double dgpatom = GetDistance(row->data[0], row->data[1], row->data[2], molecule->atoms[k].coord.x, molecule->atoms[k].coord.y, molecule->atoms[k].coord.z);
          dgpatom += pradius; /* Add the probe radius */
          dgpatom += pdistance; /* Add a probe distance if is not in direct contact */
          /*Non bondend interactions */
          /* E_qq = q1*q2/ (4 pi Epsilon0 * r) */
          double pcharge = ff.params[probe_id].charge;
          double acharge = molecule->atoms[k].charge /*ff.params[molecule->atoms[k].ffatompos].charge*/;
          el_pot += au2C(pcharge)*au2C(acharge)/(4*_pi_*epsilon0_C_J_m*dgpatom*1E-10);  // Joule
          /* E_LJ = 4 Epsilon0 * (atom_radius/r)^12 - (atom_radius/r)^6 */
          /*A and B calculated according the Lorentz-Berthelot mixing rules*/
          double r_eq_i = ff.params[molecule->atoms[k].ffatompos].vdw_R_eq;
          double r_eq_j = ff.params[probe_id].vdw_R_eq;

          double r_eq_ij = 0.5*(r_eq_i + r_eq_j);

          double eps_i = ff.params[molecule->atoms[k].ffatompos].vdw_Epsilon;
          double eps_j = ff.params[probe_id].vdw_Epsilon;

          double eps_ij = sqrt(eps_i*eps_j);

          double A = eps_ij*pow((r_eq_ij), 12);
          double B = 2.f*eps_ij*pow((r_eq_ij), 6);
          //lj_pot += 4*epsilon0_kJ_Mol*((A/pow(dgpatom, 12))-(B/pow(dgpatom, 6))); // kJoule
          lj_pot += 4*epsilon0_C_J_m*((A/pow(dgpatom, 12))-(B/pow(dgpatom, 6))); // Joule
        }

        el_pot /= 6.947857854533377e-21; // Conversion from J to kcal/mol  4184 J / 6.022E23 mol (avogadro constant)
        //lj_pot /= 4.184; // Conversion  to kcal/mol 1 kcal = 4.184 kJ
        lj_pot /= 6.947857854533377e-21; // Conversion from J to kcal/mol  4184 J / 6.022E23 mol (avogadro constant)

        /* Hydrogen bond potential calcualted according to Israelachveli book.
         * Hydrogen bond is an electrostatic (Coulombic) interaction.
         * w(r) = -Q(H+)*u * cos(t)/(4*pi*E0*E*r^2)
         * where:
         * Q(H+) is the hydrogen atom charge
         * u is the bond dipole moment
         * t is the angle formed between the dipole and the H atom
         *
         * u is calculated as follow:
         * u = pc(+)-pc(-)*dst
         * where:
         *    pc(+) is the partial charge of the first atom
         *    pc(-) is the partial charge of the second  atom
         */
        double hb_pot = 0.f;

        /*Search for hydrogen bond acceptors and donors interactions */
        if(HBACheck(patomtype) == 0){ /*Probe must be an hydrogen bond acceptor */
           POINT probecoord;
           probecoord.x = row->data[0]; probecoord.y = row->data[1]; probecoord.z = row->data[2];

           /* Interaction type
            * Probe: HBA (charge)
            * Molecule: M-OH, M-NH (dipole)
            */
          for(k = 0; k < molecule->n_bonds; k++){
            int o_id = molecule->bonds[k].origin_atom_id;
            int t_id = molecule->bonds[k].target_atom_id;
            /*char *atomsym_orig = molecule->atoms[o_id].asymbl;
            char *atomsym_targ = molecule->atoms[t_id].asymbl;*/
            if(o_id == j){  /* WE CALCULATE THE HYDROGEN BOND FOR THE (CURRENT) ATOM DEFINED AT LINE 131 */
              if(bondHBDGroup((*molecule), o_id, t_id) == 1
                /*(strcmp(atomsym_orig, "H") == 0 && strcmp(atomsym_targ, "O") == 0) ||
                 (strcmp(atomsym_orig, "H") == 0 && strcmp(atomsym_targ, "N") == 0)*/){
                /* Dipole found in molecule! */
                double u = ff.params[molecule->atoms[t_id].ffatompos].hb_Umoment;
                /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
                hb_pot += HBEnergyCalc(molecule->atoms[t_id].coord, molecule->atoms[o_id].coord, u, probecoord, ff.params[probe_id].charge);
              }
              else if(bondHBDGroup((*molecule), t_id, o_id) == 1/*(strcmp(atomsym_orig, "O") == 0 && strcmp(atomsym_targ, "H") == 0) ||
                      (strcmp(atomsym_orig, "N") == 0 && strcmp(atomsym_targ, "H") == 0)*/){ /* Atom hydrogen bond donors in molecule */
                /* Dipole found in molecule! */
                double u = ff.params[molecule->atoms[o_id].ffatompos].hb_Umoment;
                /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
                hb_pot += HBEnergyCalc(molecule->atoms[o_id].coord, molecule->atoms[t_id].coord, u, probecoord, ff.params[probe_id].charge);
              }
            }
            else if(t_id == j){
              if(bondHBDGroup((*molecule), o_id, t_id) == 1/*(strcmp(atomsym_orig, "O") == 0 && strcmp(atomsym_targ, "H")) ||
                 (strcmp(atomsym_orig, "N") == 0 && strcmp(atomsym_targ, "H"))*/){ /* Atom hydrogen bond donors in molecule */
                /* Dipole found in molecule! */
                double u = ff.params[molecule->atoms[t_id].ffatompos].hb_Umoment;
                /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
                hb_pot += HBEnergyCalc(molecule->atoms[t_id].coord, molecule->atoms[o_id].coord, u, probecoord, ff.params[probe_id].charge);
              }
              else if(bondHBDGroup((*molecule), t_id, o_id) == 1/*(strcmp(atomsym_targ, "O") == 0 && strcmp(atomsym_orig, "H")) ||
                      (strcmp(atomsym_targ, "N") == 0 && strcmp(atomsym_orig, "H"))*/){
                /* Dipole found in molecule! */
                double u = ff.params[molecule->atoms[o_id].ffatompos].hb_Umoment;
                /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
                hb_pot += HBEnergyCalc(molecule->atoms[o_id].coord, molecule->atoms[t_id].coord, u, probecoord, ff.params[probe_id].charge);
              }
            }
          }

          /* Interaction type
           * Probe: HBD (dipole)
           * Molecule: M=O, M-N, M-S (charge)
           */

          if(atomHBAGroup(patomtype) == 0){ /* Check if the probe is an hydrogen bond donor */
            for(k = 0; k < molecule->n_bonds; k++){
              size_t o_id = molecule->bonds[k].origin_atom_id;
              size_t t_id = molecule->bonds[k].target_atom_id;
              char *atomsym_orig = molecule->atoms[o_id].asymbl;
              char *atomsym_targ = molecule->atoms[t_id].asymbl;
              if(o_id == j){  /* WE CALCULATE THE HYDROGEN BOND FOR THE (CURRENT) ATOM DEFINED AT LINE 131 */
                if(atomHBAGroup(atomsym_orig) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
                  /* Dipole found in probe! */
                  double u = ff.params[probe_id].hb_Umoment;
                  double mol_charge = ff.params[molecule->atoms[o_id].ffatompos].charge;
                  /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
                  hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[o_id].coord, mol_charge);
                }

                if(atomHBAGroup(atomsym_targ) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
                  double u = ff.params[probe_id].hb_Umoment;
                  double mol_charge = ff.params[molecule->atoms[t_id].ffatompos].charge;
                  /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
                  hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[t_id].coord, mol_charge);
                }
              }
              else if(t_id == j){
                if(atomHBAGroup(molecule->atoms[o_id].type) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
                  /* Dipole found in probe! */
                  double u = ff.params[probe_id].hb_Umoment;
                  double mol_charge = ff.params[molecule->atoms[o_id].ffatompos].charge;
                  /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
                  hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[o_id].coord, mol_charge);
                }

                if(atomHBAGroup(atomsym_targ) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
                  double u = ff.params[probe_id].hb_Umoment;
                  double mol_charge = ff.params[molecule->atoms[t_id].ffatompos].charge;
                  /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
                  hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[t_id].coord, mol_charge);
                }
              }
            }
          }
        }

        hb_pot /= 6.947857854533377e-21; //conversion from Joule to kcal/mol
        //printf("%e %f\n", el_pot, ljpot);
        row->data[3] = el_pot+lj_pot+hb_pot;
        MatrixAppendRow(field, row);
        /*i++; Generate random points while increment */
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
  //DelMatrix(&angles);
  free(patomtype);
}

int PointIsCompenetrating(MOLECULE molecule, double px, double py, double pz)
{
  size_t i;
  for(i = 0; i < molecule.n_atoms; i++){
    double dp_atom = GetDistance(px, py, pz,
                                 molecule.atoms[i].coord.x,
                                 molecule.atoms[i].coord.y,
                                 molecule.atoms[i].coord.z);
    if(dp_atom < molecule.atoms[i].radius){
    /*this point is penetrating an other atom so is not in the surface*/
      return 1;
    }
    else{
     continue;
    }
  }
  return 0;
}
double Energy(MOLECULE *molecule, double px, double py, double pz, int probe_id, ForceField ff, enum RADIUS_TYPE rtype)
{
  size_t i, j;
  double pradius;
  char *patomtype = strdup(ff.params[probe_id].atom_name);
  if(rtype == vanderwaals){
    pradius = getVanDerWaalsRadiifromAtomName(patomtype);
  }
  else{
    pradius = getCovRadiusfromAtomName(patomtype);
  }

  /*Check if the voxel point is compenetrating the molecular volume */
  if(PointIsCompenetrating((*molecule), px, py, pz) == 0){ /*not compenetrating*/
    /* This is a right voxel point
      * We can calculate the potential of this point function of other atoms
      * cotributions!
      *
      * N.B.: Charges must be in atomic unit (au)!
      * 1 hartree (Eh) = 27.2113834 eV
      * 1 Hartree = 627.503 kcal/mol
      * 1 eV = 23.0605 kcal/mol
      * 1 eV = 1 J/C
      * 1 Debye = 3.336·10^(-30) C·m
      */
    /*electrostatic energy due to the charge-charge interaction*/
    double el_pot = 0.f;
    /* Lenard-Jones potential calcualted according to the classical equation*/
    double lj_pot = 0.f;
    /* Hydrogen bond potential */
    double hb_pot = 0.f;
    /* When calculate the energy of interaction we take into account
    * all the atoms
    */
    for(i = 0; i < molecule->n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dgpatom = GetDistance(px, py, pz,
                                    molecule->atoms[i].coord.x,
                                    molecule->atoms[i].coord.y,
                                    molecule->atoms[i].coord.z);
      dgpatom -= pradius; /* Remove the probe radius */
      dgpatom -= molecule->atoms[i].radius; /* Remove the atom radius */
      /*Non bondend interactions */
      /* E_qq = q1*q2/ (4 pi Epsilon0 * r) */
      double pcharge = ff.params[probe_id].charge;
      double acharge = molecule->atoms[i].charge;
      el_pot += (au2C(pcharge)*au2C(acharge))/(4*_pi_*epsilon0_C_J_m*dgpatom*1E-10); /* Joule */
      /* E_LJ = 4 Epsilon0 * (atom_radius/r)^12 - (atom_radius/r)^6 */
      /*A and B calculated according the Lorentz-Berthelot mixing rules*/
      double r_eq_i = ff.params[molecule->atoms[i].ffatompos].vdw_R_eq;
      double r_eq_j = ff.params[probe_id].vdw_R_eq;
      double r_eq_ij = 0.5*(r_eq_i + r_eq_j);
      double eps_i = ff.params[molecule->atoms[i].ffatompos].vdw_Epsilon;
      double eps_j = ff.params[probe_id].vdw_Epsilon;
      double eps_ij = sqrt(eps_i*eps_j);
      double A = eps_ij*pow((r_eq_ij), 12);
      double B = 2.f*eps_ij*pow((r_eq_ij), 6);
      /*lj_pot += 4*epsilon0_kJ_Mol*((A/pow(dgpatom,12))-(B/pow(dgpatom,6))); kJoule*/
      lj_pot += 4*epsilon0_C_J_m*((A/pow(dgpatom, 12))-(B/pow(dgpatom, 6))); /* Joule */

    /* Hydrogen bond potential calcualted according a modified version of Israelachveli book.
      * Hydrogen bond is an electrostatic (Coulombic) interaction betweem a dipole and a charge + a switch function.
      * w(r) = -Q(H+)*u * cos(t)/(4*pi*E0*E*r^2)  *switch(distance,angle)
      * where:
      * Q(H+) is the hydrogen atom charge
      * u is the bond dipole moment
      * t is the angle formed between the dipole and the H atom
      *
      * u is calculated as follow:
      * u = pc(+)-pc(-)*dst
      * or could be estimated by forcefield
      * where:
      *    pc(+) is the partial charge of the first atom
      *    pc(-) is the partial charge of the second  atom
      */

      POINT probecoord;
      probecoord.x = px; probecoord.y = py; probecoord.z = pz;
      /*Search for hydrogen bond acceptors and donors interactions */
      if(HBACheck(patomtype) == 0){ /* Probe must be an hydrogen bond acceptor */
        /* Interaction type
        * Probe: HBA (charge)
        * Molecule: M-OH, M-NH (dipole)
        *
        * Dipole considered:
        * N-H
        * O-H
        */
        for(j = 0; j < molecule->n_bonds; j++){
          int o_id = molecule->bonds[j].origin_atom_id;
          int t_id = molecule->bonds[j].target_atom_id;
          char *atomsym_orig = molecule->atoms[o_id].asymbl;
          char *atomsym_targ = molecule->atoms[t_id].asymbl;
          if(o_id == i){  /* WE CALCULATE THE HYDROGEN BOND FOR THE (CURRENT) ATOM i */
            if((strcmp(atomsym_orig, "H") == 0 && strcmp(atomsym_targ, "O") == 0) ||
                (strcmp(atomsym_orig, "H") == 0 && strcmp(atomsym_targ, "N") == 0)){
              /* Dipole found in molecule! */
              double u = ff.params[molecule->atoms[t_id].ffatompos].hb_Umoment;
              /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, d vector *probecoord, double pcharge)*/
              hb_pot += HBEnergyCalc(molecule->atoms[t_id].coord, molecule->atoms[o_id].coord, u, probecoord, ff.params[probe_id].charge);
            }
            else if((strcmp(atomsym_orig, "O") == 0 && strcmp(atomsym_targ, "H") == 0) ||
                    (strcmp(atomsym_orig, "N") == 0 && strcmp(atomsym_targ, "H") == 0)){
              /* Atom hydrogen bond donors in molecule */
              /* Dipole found in molecule! */
              double u = ff.params[molecule->atoms[o_id].ffatompos].hb_Umoment;
              /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
              hb_pot += HBEnergyCalc(molecule->atoms[o_id].coord, molecule->atoms[t_id].coord, u, probecoord, ff.params[probe_id].charge);
            }
          }
          else if(t_id == i){
            if((strcmp(atomsym_orig, "O") == 0 && strcmp(atomsym_targ, "H"))   ||
                (strcmp(atomsym_orig, "N") == 0 && strcmp(atomsym_targ, "H"))){
              /* Atom hydrogen bond donors in molecule */
              /* Dipole found in molecule! */
              double u = ff.params[molecule->atoms[o_id].ffatompos].hb_Umoment;
              /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
              hb_pot += HBEnergyCalc(molecule->atoms[o_id].coord, molecule->atoms[t_id].coord, u, probecoord, ff.params[probe_id].charge);
            }
            else if((strcmp(atomsym_targ, "O") == 0 && strcmp(atomsym_orig, "H")) ||
                    (strcmp(atomsym_targ, "N") == 0 && strcmp(atomsym_orig, "H"))){
              /* Dipole found in molecule! */
              double u = ff.params[molecule->atoms[t_id].ffatompos].hb_Umoment;
              /*double HBEnergyCalc(POINT hetatm, POINT h_atom, double u, dvector *probecoord, double pcharge)*/
              hb_pot += HBEnergyCalc(molecule->atoms[t_id].coord, molecule->atoms[o_id].coord, u, probecoord, ff.params[probe_id].charge);
            }
          }
        }
      }

    /* Interaction type
      * Probe: HBD (dipole)
      * Molecule: M=O, M-N, M-S (charge)
      */
    /*Search for hydrogen bond donors interactions */
      if(HBDCheck(patomtype) == 0){ /* Check if the probe is an hydrogen bond donor */
        for(j = 0; j < molecule->n_bonds; j++){
          size_t o_id = molecule->bonds[j].origin_atom_id;
          size_t t_id = molecule->bonds[j].target_atom_id;
          char *atomsym_orig = molecule->atoms[o_id].asymbl;
          char *atomsym_targ = molecule->atoms[t_id].asymbl;
          if(o_id == i){  /* WE CALCULATE THE HYDROGEN BOND FOR THE (CURRENT) ATOM */
            if(HBACheck(atomsym_orig) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
              /* Dipole found in probe! */
              double u = ff.params[probe_id].hb_Umoment;
              double mol_charge = ff.params[molecule->atoms[o_id].ffatompos].charge;
              /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
              hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[o_id].coord, mol_charge);
            }
            else if(HBACheck(atomsym_targ) == 0){ /* Check if the molecule is an   hydrogen bond acceptor */
              double u = ff.params[probe_id].hb_Umoment;
              double mol_charge = ff.params[molecule->atoms[t_id].ffatompos].charge;
              /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
              hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[t_id].coord, mol_charge);
            }
            else{
              continue;
            }
          }
          else if(t_id == i){
            if(HBACheck(molecule->atoms[o_id].type) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
              /* Dipole found in probe! */
              double u = ff.params[probe_id].hb_Umoment;
              double mol_charge = ff.params[molecule->atoms[o_id].ffatompos].charge;
              /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
              hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[o_id].coord, mol_charge);
            }
            else if(HBACheck(atomsym_targ) == 0){ /* Check if the molecule is an hydrogen bond acceptor */
              double u = ff.params[probe_id].hb_Umoment;
              double mol_charge = ff.params[molecule->atoms[t_id].ffatompos].charge;
              /* double HBAEnergyCalc(POINT probecoord, double u, POINT mol_coord, double mol_charge) */
              hb_pot += HBAEnergyCalc(patomtype, probecoord, u, molecule->atoms[t_id].coord, mol_charge);
            }
          }
        }
      }
    }

    /*Convert all the potential to the same scale */
    /* Conversion from J to kcal/mol  4184 J / 6.022E23 mol (avogadro constant)*/
    el_pot /= 6.947857854533377e-21;
    /* lj_pot /= 4.184;  Conversion  to kcal/mol 1 kcal = 4.184 kJ */
    /* Conversion from J to kcal/mol  4184 J / 6.022E23 mol (avogadro constant) */
    lj_pot /= 6.947857854533377e-21;
    hb_pot /= 6.947857854533377e-21; /* conversion from Joule to kcal/mol */
    free(patomtype);
    return el_pot + lj_pot + hb_pot;
  }
  else{
    free(patomtype);
    return 0.f;
  }
}

void VoxelFieldCalculator(MOLECULE *molecule, int formal_charge, int npnt, int voxel_size, enum RADIUS_TYPE rtype, int probe_id, ForceField ff, VOXEL **field)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  AssignParams(molecule, ff, rtype, formal_charge);
  //printf("Sum of charges: %.8f\n", GetTotalCharge((*molecule)));

  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;
  NewVoxel(field, npnt, npnt, npnt);
  (*field)->xmin = xmin;
  (*field)->ymin = ymin;
  (*field)->zmin = zmin;
  (*field)->xmax = xmax;
  (*field)->ymax = ymax;
  (*field)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*field)->nx-1);
  dy = (ymax-ymin)/((*field)->ny-1);
  dz = (zmax-zmin)/((*field)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*field)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*field)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*field)->nz; k++){
        z = zmin + k*dz;
        double e =  Energy(molecule, x, y, z, probe_id, ff, rtype);
        (*field)->pnt[i][j][k] = e;
      }
    }
  }
}

double ElectrostaticPotential(MOLECULE molecule, double px, double py, double pz)
{
  /*
   * E = kq/r
   * where k is the Coulomb constant
   * EPot = Sum E of each atom.
   */
  int i;
  double el_pot = 0.f;

  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    for(i = 0; i < molecule.n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      dp_atom *= 1E-10; /*Convert the distance from angstrom to meter*/
      /* calculate the potential in Joules */
      el_pot += (epsilon0_C_J_m*au2C(molecule.atoms[i].charge))/dp_atom;
    }
    /* Conversion from J to kcal/mol dividing the value for
    * 4184 J / 6.022E23 mol (avogadro constant)
    */
    el_pot /= 6.947857854533377e-21;
    return el_pot;
  }
  else{
    return 0.f;
  }
}


double InteractionPotential(MOLECULE molecule, double charge, double px, double py, double pz)
{
  /*
   * E = kq/r
   * where k is the Coulomb constant
   * EPot = Sum E of each atom.
   */
  int i;
  double el_pot = 0.f;

  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    for(i = 0; i < molecule.n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      dp_atom *= 1E-10; /*Convert the distance from angstrom to meter*/
      /* calculate the potential in Joules */
      // el_pot += (au2C(charge)*au2C(molecule.atoms[i].charge))/dp_atom;
      el_pot += (au2C(charge)*au2C(molecule.atoms[i].charge))/(4*_pi_*epsilon0_C_J_m*dp_atom); /* Joule */
    }
    /* Conversion from J to kcal/mol dividing the value for
    * 4184 J / 6.022E23 mol (avogadro constant)
    */
    el_pot /= 6.947857854533377e-21;
    return el_pot;
  }
  else{
    return 0.f;
  }
}

void VoxelElectrostaticPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **epot)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  //printf("Sum of charges: %.8f\n", GetTotalCharge((*molecule)));

  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(epot, npnt, npnt, npnt);
  (*epot)->xmin = xmin;
  (*epot)->ymin = ymin;
  (*epot)->zmin = zmin;
  (*epot)->xmax = xmax;
  (*epot)->ymax = ymax;
  (*epot)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*epot)->nx-1);
  dy = (ymax-ymin)/((*epot)->ny-1);
  dz = (zmax-zmin)/((*epot)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*epot)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*epot)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*epot)->nz; k++){
        z = zmin + k*dz;
        (*epot)->pnt[i][j][k] = ElectrostaticPotential((*molecule), x, y, z);
      }
    }
  }
}

void VoxelElectrostaticPotentialInteractionCalculator(MOLECULE *molecule,
                                                      double charge,
                                                      int npnt,
                                                      int voxel_size,
                                                      enum RADIUS_TYPE rtype,
                                                      VOXEL **epot)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  //printf("Sum of charges: %.8f\n", GetTotalCharge((*molecule)));

  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(epot, npnt, npnt, npnt);
  (*epot)->xmin = xmin;
  (*epot)->ymin = ymin;
  (*epot)->zmin = zmin;
  (*epot)->xmax = xmax;
  (*epot)->ymax = ymax;
  (*epot)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*epot)->nx-1);
  dy = (ymax-ymin)/((*epot)->ny-1);
  dz = (zmax-zmin)/((*epot)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*epot)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*epot)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*epot)->nz; k++){
        z = zmin + k*dz;
        (*epot)->pnt[i][j][k] = InteractionPotential((*molecule), charge, x, y, z);
      }
    }
  }
}

void ProteinLigandVoxelElectrostaticPotentialCalculator(MOLECULE *protein,
                                                        MOLECULE *ligand,
                                                        int npnt,
                                                        int voxel_size,
                                                        enum RADIUS_TYPE rtype,
                                                        VOXEL **pocket_epot,
                                                        VOXEL **ligand_epot)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  //printf("Sum of charges: %.8f\n", GetTotalCharge((*molecule)));

  // Get the center based on the ligand
  GetMolBox((*ligand), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(pocket_epot, npnt, npnt, npnt);
  (*pocket_epot)->xmin = xmin;
  (*pocket_epot)->ymin = ymin;
  (*pocket_epot)->zmin = zmin;
  (*pocket_epot)->xmax = xmax;
  (*pocket_epot)->ymax = ymax;
  (*pocket_epot)->zmax = zmax;

  NewVoxel(ligand_epot, npnt, npnt, npnt);
  (*ligand_epot)->xmin = xmin;
  (*ligand_epot)->ymin = ymin;
  (*ligand_epot)->zmin = zmin;
  (*ligand_epot)->xmax = xmax;
  (*ligand_epot)->ymax = ymax;
  (*ligand_epot)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < ligand->n_atoms; i++){
    if(rtype == vanderwaals){
      ligand->atoms[i].radius = getVanDerWaalsRadiifromAtomName(ligand->atoms[i].asymbl);
    }
    else{
      ligand->atoms[i].radius = getCovRadiusfromAtomName(ligand->atoms[i].asymbl);
    }
  }

  for(i = 0; i < protein->n_atoms; i++){
    if(rtype == vanderwaals){
      protein->atoms[i].radius = getVanDerWaalsRadiifromAtomName(protein->atoms[i].asymbl);
    }
    else{
      protein->atoms[i].radius = getCovRadiusfromAtomName(protein->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*ligand_epot)->nx-1);
  dy = (ymax-ymin)/((*ligand_epot)->ny-1);
  dz = (zmax-zmin)/((*ligand_epot)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*ligand_epot)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*ligand_epot)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*ligand_epot)->nz; k++){
        z = zmin + k*dz;
        (*pocket_epot)->pnt[i][j][k] = ElectrostaticPotential((*protein), x, y, z);
        (*ligand_epot)->pnt[i][j][k] = ElectrostaticPotential((*ligand), x, y, z);
      }
    }
  }
}

double VDWPotential(MOLECULE molecule, double px, double py, double pz)
{
  int i;
  double e_vdw = 0.f;

  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    for(i = 0; i < molecule.n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      /* dp_atom is already in angstrom to meter*/
      /* calculate the potential in kcal/mol */
      double D_ii, x_ii;
      getVanDerWaalsParams(molecule.atoms[i].asymbl, &D_ii, &x_ii);
      /*
       * See the UFF VdW potential in doi:10.1.1.208.7677
       * UFF, a Full Periodic Table Force Field for
       * Molecular Mechanics and Molecular Dynamics Simulations
       */
       e_vdw += D_ii*(-2.f*pow(x_ii/dp_atom, 6)+pow(x_ii/dp_atom, 12));
    }
    return e_vdw;
  }
  else{
    return 0.f;
  }
}

void VoxelVdWPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **vdwp)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(vdwp, npnt, npnt, npnt);
  (*vdwp)->xmin = xmin;
  (*vdwp)->ymin = ymin;
  (*vdwp)->zmin = zmin;
  (*vdwp)->xmax = xmax;
  (*vdwp)->ymax = ymax;
  (*vdwp)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*vdwp)->nx-1);
  dy = (ymax-ymin)/((*vdwp)->ny-1);
  dz = (zmax-zmin)/((*vdwp)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*vdwp)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*vdwp)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*vdwp)->nz; k++){
        z = zmin + k*dz;
        (*vdwp)->pnt[i][j][k] = VDWPotential((*molecule), x, y, z);
      }
    }
  }
}

double ENegPotential(MOLECULE molecule, double px, double py, double pz)
{
  /*
   * Like in MLP the distance function
   * is written as follow:
   * EnegP_k = Sum(atom from i to N) Eneg_i * exp(-d_ik/2)
   * Where:
   * d_ik is the distance of atom i to the point k
   *   (Fauchere doi:10.1016/S0263-7855(98)80004-0)
   * Eneg_i is the eletronegativity of the i atom
   */
  int i;
  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    double e_neg_p = 0.f;
    for(i = 0; i < molecule.n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      double eneg_i = getEnegativity(molecule.atoms[i].asymbl);
      e_neg_p += eneg_i * exp(-1*dp_atom/2.f);
    }
    return e_neg_p;
  }
  else{
    return 0.f;
  }
}

void VoxelENegPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **enegp)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(enegp, npnt, npnt, npnt);
  (*enegp)->xmin = xmin;
  (*enegp)->ymin = ymin;
  (*enegp)->zmin = zmin;
  (*enegp)->xmax = xmax;
  (*enegp)->ymax = ymax;
  (*enegp)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*enegp)->nx-1);
  dy = (ymax-ymin)/((*enegp)->ny-1);
  dz = (zmax-zmin)/((*enegp)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*enegp)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*enegp)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*enegp)->nz; k++){
        z = zmin + k*dz;
        (*enegp)->pnt[i][j][k] = ENegPotential((*molecule), x, y, z);
      }
    }
  }
}


double LJPotential(MOLECULE molecule,
                        AtomsProperty *epsilon,
                        AtomsProperty *sigma,
                        double px, double py, double pz)
{
  /*
   * Generate the Lennard-Jones potential epsilon and sigma properties
   * V(r) = 4epsilon * [(sigma/r)^12 - (sigma/r)^6]
   * Where:
   * r is the distance of atom i to the point k
   * epsilon is the potential depth
   * sigma is the sphere diameter
   * ConstAtom_i is the empiric constant for the i atom
   *
   */
  int i;
  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    double pot = 0.f;
    for(i = 0; i < molecule.n_atoms; i++){
      /* we add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      double eps_i = getGenericAProperty(molecule.atoms[i].asymbl, epsilon);
      double sigma_i = getGenericAProperty(molecule.atoms[i].asymbl, sigma);

      /*
      * Normally to reduce the computational time, the L-J interaction
      * is computed when the distance r divided with the sphere diameter sigma
      * is <= 2.5. This approximation is known to be the "Truncated (and shifted) form
      * of L-J potential"
      *
      * if(dp_atom/sigma_i <= 2.5){
      *   pot += 4.f*eps_i* (pow((sigma_i/dp_atom), 12) - pow((sigma_i/dp_atom), 6));
      * }
      * else{
      *     return pot;
      * }
      */
      pot += 4.f*eps_i* (pow((sigma_i/dp_atom), 12) - pow((sigma_i/dp_atom), 6));
    }
    return pot;
  }
  else{
    return 0.f;
  }
}

void VoxelLJPotentialCalculator(MOLECULE *molecule,
                                AtomsProperty *epsilon,
                                AtomsProperty *sigma,
                                int npnt,
                                int voxel_size,
                                enum RADIUS_TYPE rtype,
                                VOXEL **ljpot)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(ljpot, npnt, npnt, npnt);
  (*ljpot)->xmin = xmin;
  (*ljpot)->ymin = ymin;
  (*ljpot)->zmin = zmin;
  (*ljpot)->xmax = xmax;
  (*ljpot)->ymax = ymax;
  (*ljpot)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*ljpot)->nx-1);
  dy = (ymax-ymin)/((*ljpot)->ny-1);
  dz = (zmax-zmin)/((*ljpot)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*ljpot)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*ljpot)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*ljpot)->nz; k++){
        z = zmin + k*dz;
        (*ljpot)->pnt[i][j][k] = LJPotential((*molecule), epsilon, sigma, x, y, z);
      }
    }
  }
}

double GenericPotential(MOLECULE molecule, AtomsProperty *lst, double px, double py, double pz)
{
  /*
   * Like in MLP the distance function
   * is written as follow:
   * GenP_k = Sum(atom from i to N) ConstAtom_i * exp(-d_ik/2)
   * Where:
   * d_ik is the distance of atom i to the point k
   *   (Fauchere doi:10.1016/S0263-7855(98)80004-0)
   * ConstAtom_i is the empiric constant for the i atom
   */
  int i;
  /*If the point is not compenetrating any atoms*/
  if(PointIsCompenetrating(molecule, px, py, pz) == 0){
    double pot = 0.f;
    for(i = 0; i < molecule.n_atoms; i++){
      /* we should add the atom probe radius to simulate the sphere contact*/
      double dp_atom = GetDistance(px, py, pz,
                                   molecule.atoms[i].coord.x,
                                   molecule.atoms[i].coord.y,
                                   molecule.atoms[i].coord.z);
      double pot_i = getGenericAProperty(molecule.atoms[i].type, lst);
      pot += pot_i * exp(-1*dp_atom/2.f);
      //pot += pot_i/dp_atom;
    }
    return pot;
  }
  else{
    return 0.f;
  }
}

void VoxelGenericPotentialCalculator(MOLECULE *molecule, AtomsProperty *lst,  int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **gpot)
{
  int i, j, k;
  double  xmin, xmax, ymin, ymax, zmin, zmax;
  double x, y, z, dx, dy, dz;

  /*AssignParams(molecule, ff, rtype, formal_charge);*/
  GetMolBox((*molecule), (double)voxel_size, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

  /*
  double ang2au = 1.0 / 0.5291772083; //0.5291772083 is the bohr radius
  xmin *= ang2au;
  xmax *= ang2au;
  ymin *= ang2au;
  ymax *= ang2au;
  zmin *= ang2au;
  zmax *= ang2au;*/

  //printf("xmin %f xmax %f ymin %f ymax %f zmin %f zmax %f\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /*Create Voxel*/
  dx = (xmax-xmin)/(double)npnt;
  dy = (ymax-ymin)/(double)npnt;
  dz = (zmax-zmin)/(double)npnt;

  NewVoxel(gpot, npnt, npnt, npnt);
  (*gpot)->xmin = xmin;
  (*gpot)->ymin = ymin;
  (*gpot)->zmin = zmin;
  (*gpot)->xmax = xmax;
  (*gpot)->ymax = ymax;
  (*gpot)->zmax = zmax;

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /* Calculate the field in each point of the voxel */

  dx = (xmax-xmin)/((*gpot)->nx-1);
  dy = (ymax-ymin)/((*gpot)->ny-1);
  dz = (zmax-zmin)/((*gpot)->nz-1);

  /* create voxel point using an equally-spaced numbers algoritm
   * division first: start + i*(stop-start)/(num-1)
   */
  for(i = 0; i < (*gpot)->nx; i++){
    x =  xmin + i*dx;
    for(j = 0; j < (*gpot)->ny; j++){
      y = ymin + j*dy;
      for(k = 0; k < (*gpot)->nz; k++){
        z = zmin + k*dz;
        (*gpot)->pnt[i][j][k] = GenericPotential((*molecule), lst, x, y, z);
      }
    }
  }
}


/*
 * Calculate the electrostatic potential using the point charge approach
 * and the partial charges of the molecules. The Electrostatic Potential
 * is expressed in kcal/mol
 */
void SphericalElectrostaticPotentialCalculator(MOLECULE *molecule,
                                               int npnt,
                                               enum RADIUS_TYPE rtype,
                                               matrix *epot){
  size_t i, j;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  matrix *pnts;

  NewDVector(&row, 4);

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);
    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */

      row->data[0] = pnts->data[i][0] + molecule->atoms[j].coord.x;
      row->data[1] = pnts->data[i][1] + molecule->atoms[j].coord.y;
      row->data[2] = pnts->data[i][2] + molecule->atoms[j].coord.z;

      if(PointIsCompenetrating((*molecule), row->data[0], row->data[1], row->data[2]) == 0){ /*not compenetrating*/
        row->data[3] = ElectrostaticPotential((*molecule), row->data[0], row->data[1], row->data[2]);
        MatrixAppendRow(epot, row);
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
}

/*
 * Calculate the Van der Waals potential using the point charge approach
 * and the UFF VdW parameters.
 * The potential is expressed in kcal/mol
 */
void SphericalVdWPotentialCalculator(MOLECULE *molecule,
                                     int npnt,
                                     enum RADIUS_TYPE rtype,
                                     matrix *vdwp){
  size_t i, j;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  matrix *pnts;

  NewDVector(&row, 4);

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);
    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */

      row->data[0] = pnts->data[i][0] + molecule->atoms[j].coord.x;
      row->data[1] = pnts->data[i][1] + molecule->atoms[j].coord.y;
      row->data[2] = pnts->data[i][2] + molecule->atoms[j].coord.z;

      if(PointIsCompenetrating((*molecule), row->data[0], row->data[1], row->data[2]) == 0){ /*not compenetrating*/
        row->data[3] = VDWPotential((*molecule), row->data[0], row->data[1], row->data[2]);
        MatrixAppendRow(vdwp, row);
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
}

/*
 * Calculate the Electronegativity potential (distribution)
 * using the point charge approach
 * The potential is expressed in pauling relative electronegativity scale
 * It seems to be correlated with the molecular shape...
 */
void SphericalENegPotentialCalculator(MOLECULE *molecule, int npnt, enum RADIUS_TYPE rtype, matrix *enegp){
  size_t i, j;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  matrix *pnts;

  NewDVector(&row, 4);

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);
    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */

      row->data[0] = pnts->data[i][0] + molecule->atoms[j].coord.x;
      row->data[1] = pnts->data[i][1] + molecule->atoms[j].coord.y;
      row->data[2] = pnts->data[i][2] + molecule->atoms[j].coord.z;

      if(PointIsCompenetrating((*molecule), row->data[0], row->data[1], row->data[2]) == 0){ /*not compenetrating*/
        row->data[3] = ENegPotential((*molecule), row->data[0], row->data[1], row->data[2]);
        MatrixAppendRow(enegp, row);
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
}

/*
 * Calculate the Lennard-Jones potential using the point charge approach
 * and epsilon sigma external properties
 * The potential is expressed in kcal/mol
 */
void SphericalLJPotentialCalculator(MOLECULE *molecule,
                                    AtomsProperty *epsilon,
                                    AtomsProperty *sigma,
                                    int npnt,
                                    enum RADIUS_TYPE rtype,
                                    matrix *ljpot){
  size_t i, j;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  matrix *pnts;

  NewDVector(&row, 4);

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);
    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */

      row->data[0] = pnts->data[i][0] + molecule->atoms[j].coord.x;
      row->data[1] = pnts->data[i][1] + molecule->atoms[j].coord.y;
      row->data[2] = pnts->data[i][2] + molecule->atoms[j].coord.z;

      if(PointIsCompenetrating((*molecule), row->data[0], row->data[1], row->data[2]) == 0){ /*not compenetrating*/
        row->data[3] = LJPotential((*molecule), epsilon, sigma, row->data[0], row->data[1], row->data[2]);
        MatrixAppendRow(ljpot, row);
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
}

/*
 * Calculate the Electronegativity potential (distribution)
 * using the point charge approach
 * The potential is expressed in pauling relative electronegativity scale
 * It seems to be correlated with the molecular shape...
 */
void SphericalGenericPotentialCalculator(MOLECULE *molecule,
                                         AtomsProperty *lst,
                                         int npnt,
                                         enum RADIUS_TYPE rtype,
                                         matrix *gpot){
  size_t i, j;
  srand((unsigned)npnt+molecule->n_atoms+molecule->n_bonds);
  dvector *row;
  matrix *pnts;

  NewDVector(&row, 4);

  /* set the radius type */
  for(i = 0; i < molecule->n_atoms; i++){
    if(rtype == vanderwaals){
      molecule->atoms[i].radius = getVanDerWaalsRadiifromAtomName(molecule->atoms[i].asymbl);
    }
    else{
      molecule->atoms[i].radius = getCovRadiusfromAtomName(molecule->atoms[i].asymbl);
    }
  }

  /*1) For each atom generate npnt random points part of the surface.
    2) If the point in surface is not overlapping an other atom surface then calcualte the energy
    3) Calculate the energy using the electrostacir potential, the lennar jones and hydrogen bond
    */
  for(j = 0; j < molecule->n_atoms; j++){
    initMatrix(&pnts);
    GenEquidistributedPointsInSphere(npnt, molecule->atoms[j].radius, pnts);
    for(i = 0; i < pnts->row; i++){
    /*  Method used with SphereCordaAngleGenerator()
    for(i = 0; i < angles->row; i++){

      double theta = angles->data[i][0];
      double phi =  angles->data[i][1];

      row->data[0] = (molecule->atoms[j].radius*cos(theta)*sin(phi)) + molecule->atoms[j].coord.x;
      row->data[1] = (molecule->atoms[j].radius*sin(theta)*sin(phi)) + molecule->atoms[j].coord.y;
      row->data[2] = (molecule->atoms[j].radius*cos(phi)) + molecule->atoms[j].coord.z;
      */

      row->data[0] = pnts->data[i][0] + molecule->atoms[j].coord.x;
      row->data[1] = pnts->data[i][1] + molecule->atoms[j].coord.y;
      row->data[2] = pnts->data[i][2] + molecule->atoms[j].coord.z;

      if(PointIsCompenetrating((*molecule), row->data[0], row->data[1], row->data[2]) == 0){ /*not compenetrating*/
        row->data[3] = GenericPotential((*molecule), lst, row->data[0], row->data[1], row->data[2]);
        MatrixAppendRow(gpot, row);
      }
      else{
        continue;
      }
    }
    DelMatrix(&pnts);
  }
  DelDVector(&row);
}

void SaveDXField(VOXEL *field, char *outfile)
{
  size_t i, j, k, n;
  FILE *fp;
  fp = fopen(outfile, "w");
  fprintf(fp, "object 1 class gridpositions counts %lu %lu %lu\n", field->nx, field->ny, field->nz);
  fprintf(fp, "origin %f %f %f\n", field->xmin, field->ymin, field->zmin);
  double hx = (field->xmax-field->xmin)/(double)field->nx;
  double hy = (field->ymax-field->ymin)/(double)field->ny;
  double hz = (field->zmax-field->zmin)/(double)field->nz;
  fprintf(fp, "delta %f 0.0 0.0\n", hx);
  fprintf(fp, "delta 0.0 %f 0.0\n", hy);
  fprintf(fp, "delta 0.0 0.0 %f\n", hz);
  fprintf(fp, "object 2 class gridconnections counts %lu %lu %lu\n", field->nx, field->ny, field->nz);
  fprintf(fp, "object 3 class array type float rank 0 items %lu data follows\n", field->nx*field->ny*field->nz);
  n = 0;
  for(i = 0; i < field->nx; i++){
    for(j = 0; j < field->ny; j++){
      for(k = 0; k < field->nz; k++){
        if(n == 2){
          fprintf(fp, "%.6e\n", field->pnt[i][j][k]);
          n = 0;
        }
        else{
          fprintf(fp, "%.6e ", field->pnt[i][j][k]);
          n++;
        }
      }
    }
  }

  if(n != 0)
    fprintf(fp, "\n");

  /*fprintf(fp, "attribute \"dep\" string \"positions\"\n");*/

  fprintf(fp, "object \"regular positions regular connections\" class field\n");
  fprintf(fp, "component \"positions\" value 1\n");
  fprintf(fp, "component \"connections\" value 2\n");
  fprintf(fp, "component \"data\" value 3");
  fclose(fp);
}

void SaveCubeFile(MOLECULE molecule, VOXEL *field, char *outfile)
{
  size_t i, j, k, n;
  double ang2au = 1.0 / 0.5291772083;
  FILE *fp;
  fp = fopen(outfile, "w");
  fprintf(fp, "CUBE FILE\n");
  fprintf(fp, "Molecular Interaction Fields for %s\n", molecule.molname);
  fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", molecule.n_atoms, field->xmin, field->ymin, field->zmin);
  fprintf(fp, "%5lu%12.6f%12.6f%12.6f\n", field->nx, (field->xmax - field->xmin)/field->nx, 0.f, 0.f);
  fprintf(fp, "%5lu%12.6f%12.6f%12.6f\n", field->ny, 0.f, (field->ymax - field->ymin)/field->ny, 0.f);
  fprintf(fp, "%5lu%12.6f%12.6f%12.6f\n", field->ny, 0.f, 0.f, (field->zmax - field->zmin)/field->nz);
  for(i = 0; i < molecule.n_atoms; i++){
    fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n", getAtomNumberfromAtomName(molecule.atoms[i].asymbl), 0.f, molecule.atoms[i].coord.x*ang2au, molecule.atoms[i].coord.y*ang2au, molecule.atoms[i].coord.z*ang2au);
  }

  n = 0;
  for(i = 0; i < field->nx; i++){
    for(j = 0; j < field->ny; j++){
      for(k = 0; k < field->nz; k++){
        fprintf(fp, "%14.5e", field->pnt[i][j][k]);
        n += 1;
        if(n > 5){
          fprintf(fp, "\n");
          n = 0;
        }
      }
      if(n != 0){
        fprintf(fp, "\n");
        n = 0;
      }
    }
  }
  fclose(fp);
}

void SaveMol2VoxelField(VOXEL *field, char *outfile)
{
 size_t i, j, k;
 FILE *fp;
 fp = fopen(outfile, "w");
 fprintf(fp, "@<TRIPOS>MOLECULE\n");
 fprintf(fp, "Volume Points\n");

 fprintf(fp, "%5zu   0 0 0 0\n", field->nx*field->ny*field->nz);
 fprintf(fp, "SMALL\n");
 fprintf(fp, "NO_CHARGES\n\n");
 fprintf(fp, "@<TRIPOS>ATOM\n");

 double dx = (field->xmax - field->xmin)/field->nx;
 double dy = (field->ymax - field->ymin)/field->ny;
 double dz = (field->zmax - field->zmin)/field->nz;

 double x = field->xmin;
 for(i = 0; i < field->nx; i++){
   double y = field->ymin;
   for(j = 0; j < field->ny; j++){
     double z = field->ymin;
     for(k = 0; k < field->nz; k++){
       fprintf(fp, "%7zu ****    %10.2f%10.2f%10.2f Du      1  UNL1       %7.4f\n", i+1, x, y, z, field->pnt[i][j][k]);
       z += dz;
     }
     y += dy;
   }
   x += dx;
 }
 fclose(fp);
}

void WriteMol2SphericalPoints(matrix *field, char *outfile)
{
 size_t i;
 FILE *fp;
 fp = fopen(outfile, "w");
 fprintf(fp, "@<TRIPOS>MOLECULE\n");
 fprintf(fp, "Volume Points\n");

 fprintf(fp, "%5zu   0 0 0 0\n", field->row);
 fprintf(fp, "SMALL\n");
 fprintf(fp, "NO_CHARGES\n\n");
 fprintf(fp, "@<TRIPOS>ATOM\n");
 for(i = 0; i < field->row; i++){
   fprintf(fp, "%7zu ****    %10.2f%10.2f%10.2f Du      1  UNL1       %7.4f\n", i+1, field->data[i][0], field->data[i][1], field->data[i][2], field->data[i][3]);
 }
 fclose(fp);
}
