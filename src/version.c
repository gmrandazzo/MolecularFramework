/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include "version.h"


void get_version(int *major, int *minor, int *patch)
{
  (*major) = (int)MAJOR;
  (*minor) = (int)MINOR;
  (*patch) = (int)PATCH;
}
