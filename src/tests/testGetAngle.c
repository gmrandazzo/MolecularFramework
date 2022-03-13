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
  POINT c;

  p1.x = -5.1933; p1.y = 1.6483; p1.z = -0.7141; //H1
  p2.x = -3.9465; p2.y = 2.1928; p2.z = -0.0208; //H2
  c.x = -4.9144; c.y = 2.2214; c.z = 0.0157; //O

  double angle = GetAngle(p1, p2, c);
  double grad_angle = (angle*360)/ (2*3.14);
  printf("Angle between H-O-H rad: %f;  °: %f\n", angle, grad_angle);
  return 0;
}
