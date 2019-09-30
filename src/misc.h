/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#ifndef MISC_H
#define MISC_H

/*
 * Misc Functions
 * Updated at 19 Oct 2016
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>


#define EPSILON 1e-3  /* Define your own tolerance*/
#define FLOAT_EQ(x,v, EPSILON) (((v - EPSILON) < x) && (x <( v + EPSILON)))

void rtrim(char *str);
void ltrim(char *str);
void trim(char *str);
char *Trim(char *s);
void copystr(char *src, char *dst);
char* concat(int count, ...);

/*
 * Calculate the square of a number
 */
double square(double x);

/*
 * Trim string
 */
char *StrTrim(char *s);

/* split a string using.
 * Example
 * char **tokens;
 *  char months[] = "JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC";
 * tokens = str_split(monts, ',');
 * free(tokens);
 *
 */
char** str_split(char* a_str, const char a_delim);

#endif
