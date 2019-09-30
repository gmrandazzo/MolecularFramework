/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef VERSION_H
#define VERSION_H

#include <stdio.h>
#include <stdlib.h>

#define MAJOR 1
#define MINOR 0
#define PATCH 0

void get_libname(char *lname);
void get_version(int *major, int *minor, int *patch);

#endif
