/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef MOL3DFIELDS_H
#define MOL3DFIELDS_H
#include <stdio.h>
#include <scientific.h>

#include "molecule.h"
#include "periodic_table.h"
#include "forcefield.h"

#define C_e 1.60217656535E-19
#define epsilon0_C_J_m 8.854187817620E-12 /* C2/J·m */
#define epsilon0_kJ_mol 0.997
#define epsilon0_kcal_mol 0.997

/* Convert from atomic units to Coulomb */
#define au2C(x) x*1.60217656535E-19
/* Convert dipole moment from debye to Coulomb*m*/
#define D2Cm(x) x*3.336E-30

typedef struct {
  double ***pnt;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  size_t nx;
  size_t ny;
  size_t nz;
} VOXEL;

void NewVoxel(VOXEL **v, size_t nx, size_t ny, size_t nz);
void DelVoxel(VOXEL **v);

/*Convert a voxel into a matrix with 4 colums:
 * col1: x coord
 * col2: y coord
 * col3: z coord
 * col4: voxel value
 */
void Voxel2Matrix(VOXEL *v, matrix *m);


/* Function to load probes and other parameters */
/*
 * Calculate the interaction fields using the voxel approach.
 *
 * Input parameters
 * molecule (struct): the molecular graph
 * size_x, size_y, size_z (double): the voxel size
 * gres (double): the grid spacing resolution (0.33 Angstrom, ...)
 * rtype (enum): radius type (van der waals or covalent)
 *
 * Output:
 * field (matrix strct): voxel output
 */
void FieldCalculator(ForceField ff,
                     MOLECULE *molecule,
                     int formal_charge,
                     size_t npnt,
                     enum RADIUS_TYPE rtype,
                     int probe_id,
                     double pdistance,
                     matrix *field);


/*
 * Calculate the electrostatic potential using the point charge approach
 * and the partial charges of the molecules. The Electrostatic Potential
 * is expressed in kcal/mol
 */
void SphericalElectrostaticPotentialCalculator(MOLECULE *molecule,
                                               int npnt,
                                               enum RADIUS_TYPE rtype,
                                               matrix *epot);

/*
 * Calculate the Van der Waals potential using the point charge approach
 * and the UFF VdW parameters.
 * The potential is expressed in kcal/mol
 */
void SphericalVdWPotentialCalculator(MOLECULE *molecule,
                                     int npnt,
                                     enum RADIUS_TYPE rtype,
                                     matrix *vdwp);


/*
 * Calculate the Electronegativity potential (distribution)
 * using the point charge approach
 * The potential is expressed in pauling relative electronegativity scale
 * It seems to be correlated with the molecular shape...
 */
void SphericalENegPotentialCalculator(MOLECULE *molecule,
                                      int npnt,
                                      enum RADIUS_TYPE rtype,
                                      matrix *enegp);


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
                                    matrix *gpot);

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
                                         matrix *enegp);


/*
 * Calculate the interaction fields using the voxel approach.
 *
 * Input parameters
 * molecule (struct): the molecular graph
 * size_x, size_y, size_z (double): the voxel size
 * gres (double): the grid spacing resolution (0.33 Angstrom, ...)
 * rtype (enum): radius type (van der waals or covalent)
 *
 * Output:
 * field (matrix strct): voxel output
 */
void VoxelFieldCalculator(MOLECULE *molecule, int formal_charge, int npnt, int voxel_size, enum RADIUS_TYPE rtype, int probe_id, ForceField ff, VOXEL **field);

/*
 * Calculate the electrostatic potential using the point charge approach
 * and the partial charges of the molecules. The Electrostatic Potential
 * is expressed in kcal/mol
 */
void VoxelElectrostaticPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **epot);

/*
 * Calculate the electrostatic potential using the point charge approach
 * and the partial charges of a ligand in a protein. You will get two
 * electrostatic potentials which are complementary and independent one each other
 * The Electrostatic Potential is expressed in kcal/mol
 */
void ProteinLigandVoxelElectrostaticPotentialCalculator(MOLECULE *protein,
                                                        MOLECULE *ligand,
                                                        int npnt,
                                                        int voxel_size,
                                                        enum RADIUS_TYPE rtype,
                                                        VOXEL **pocket_epot,
                                                        VOXEL **ligand_epot);

/*
 * Calculate the Van der Waals potential using the point charge approach
 * and the UFF VdW parameters.
 * The potential is expressed in kcal/mol
 */
void VoxelVdWPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **vdwp);


/*
 * Calculate the Electronegativity potential (distribution)
 * using the point charge approach
 * The potential is expressed in pauling relative electronegativity scale
 * It seems to be correlated with the molecular shape...
 */
void VoxelENegPotentialCalculator(MOLECULE *molecule, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **enegp);


/*
 * Calculate the Lennard-Jones potential using the point charge approach
 * and epsilon sigma external properties
 * The potential is expressed in kcal/mol
 */
void VoxelLJPotentialCalculator(MOLECULE *molecule, AtomsProperty *epsilon, AtomsProperty *sigma, int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **gpot);


/*
 * Calculate the Electronegativity potential (distribution)
 * using the point charge approach
 * The potential is expressed in pauling relative electronegativity scale
 * It seems to be correlated with the molecular shape...
 */

void VoxelGenericPotentialCalculator(MOLECULE *molecule, AtomsProperty *lst,  int npnt, int voxel_size, enum RADIUS_TYPE rtype, VOXEL **enegp);

void SaveDXField(VOXEL *fields, char *outfile);
void SaveCubeFile(MOLECULE molecule, VOXEL *fields, char *outfile);
void SaveMol2VoxelField(VOXEL *field, char *outfile);
void WriteMol2SphericalPoints(matrix *field, char *outfile);

#endif
