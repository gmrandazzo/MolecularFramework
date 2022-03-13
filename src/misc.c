/*
 * (c) 2012-2019 gmrandazzo@gmail.com
 * This file is part of MolecularFramework.
 * You can use,modify, and distribute it under
 * the terms of the GNU General Public Licenze, version 3.
 * See the file LICENSE for details
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

void rtrim(char *str)
{
  size_t n;
  n = strlen(str);
  while (n > 0 && isspace((unsigned char)str[n - 1])) {
    n--;
  }
  str[n] = '\0';
}

void ltrim(char *str)
{
  size_t n;
  n = 0;
  while (str[n] != '\0' && isspace((unsigned char)str[n])) {
    n++;
  }
  memmove(str, str + n, strlen(str) - n + 1);
}

/*
 * Leading/trailing whitespaces in strings working on the same string
 */
void trim(char *str)
{
  rtrim(str);
  ltrim(str);
}

/*
 * Leading/trailing whitespaces in strings returning a new string
 */
char *Trim(char *s)
{
    char *ptr;
    if (!s)
        return NULL;   // handle NULL string
    if (!*s)
        return s;      // handle empty string
    for (ptr = s + strlen(s) - 1; (ptr >= s) && isspace(*ptr); --ptr);
    ptr[1] = '\0';
    return s;
}

void copystr(char *src, char *dst){
  int j;
  int pos = 0;
  for (j = 0; j < strlen(src); j++) {
    pos += sprintf(&dst[pos], "%c", src[j]);
  }
}

/* Concatenate two string.
 * Usage:  char *str = concat(2,"a","b");
 * printf("%s\n", str);
 * free(str);
 * >>> ab
 */
char* concat(int count, ...)
{
  va_list ap;
  int i;

  // Find required length to store merged string
  int len = 1; // room for NULL
  va_start(ap, count);
  for(i=0 ; i<count ; i++)
    len += strlen(va_arg(ap, char*));
  va_end(ap);

  // Allocate memory to concat strings
  char *merged = calloc(sizeof(char),len);
  int null_pos = 0;

  // Actually concatenate strings
  va_start(ap, count);
  for(i=0 ; i<count ; i++){
    char *s = va_arg(ap, char*);
    strcpy(merged+null_pos, s);
    null_pos += strlen(s);
  }
  va_end(ap);

  return merged;
}

double square(double x)
{
  return x*x;
}

char *StrTrim(char *s)
{
  char *ptr;
  if (!s)
      return NULL;   // handle NULL string
  if (!*s)
      return s;      // handle empty string
  for (ptr = s + strlen(s) - 1; (ptr >= s) && isspace(*ptr); --ptr);
  ptr[1] = '\0';
  return s;
}

/* split a string using.
 * Example
 * char **tokens;
 * char months[] = "JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC";
 * tokens = str_split(monts, ',');
 *if (tokens)
 *   {
 *        int i;
 *        for (i = 0; *(tokens + i); i++)
 *        {
 *            printf("month=[%s]\n", *(tokens + i));
 *            free(*(tokens + i));
 *        }
 *       printf("\n");
 *       free(tokens);
 *    }
 *
 * free(tokens);
 *
 */
char** str_split(char* a_str, const char a_delim)
{
  char** result    = 0;
  size_t count     = 0;
  char* tmp        = a_str;
  char* last_comma = 0;
  char delim[2];
  delim[0] = a_delim;
  delim[1] = 0;

  /* Count how many elements will be extracted. */
  while (*tmp){
    if (a_delim == *tmp){
      count++;
      last_comma = tmp;
    }
    tmp++;
  }

  /* Add space for trailing token. */
  count += last_comma < (a_str + strlen(a_str) - 1);

  /* Add space for terminating null string so caller
      knows where the list of returned strings ends. */
  count++;

  result = malloc(sizeof(char*) * count);

  if (result){
    size_t idx  = 0;
    char* token = strtok(a_str, delim);

    while (token){
      assert(idx < count);
      *(result + idx++) = strdup(token);
      token = strtok(0, delim);
    }

    assert(idx == count - 1);
    *(result + idx) = 0;
  }
  return result;
}
