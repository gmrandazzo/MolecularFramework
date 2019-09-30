/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include "../geomdesc.h"

int main(void)
{
  POINT p1;
  POINT p2;
  POINT p3;
  POINT p4;


  p1.x = 1219.7; p1.y = 4106.1; p1.z = -7366.7;
  p2.x = 1277.1; p2.y = 4016.6; p2.z = -7447.1;
  p3.x = 1398.6; p3.y = 3944.8; p3.z = -7407.8;
  p4.x = 1501.2; p4.y = 3943.2; p4.z = -7521.3;


  double dihedral = DihedralAngle(p1, p2, p3, p4);
  double grad_dihedral = (dihedral*360)/ (2*3.14);
  if(FLOAT_EQ(dihedral, -2.3420649099786144, 1e-10))
    printf("Test \"DihedralAngle()\": OK!\n");
  printf("Angle between p1,p2,p3,p4 rad: %f;  °: %f\n", dihedral, grad_dihedral);
  //('Dihedral angle in radians:', -2.3420649099786144, 'and degrees:', -134.19043468746167)
  return 0;
}
