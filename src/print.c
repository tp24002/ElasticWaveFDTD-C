#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/print.h"
#include "../header/struct.h"

void print(double **f, FILE *fp, Coord Txx) {
  int i, k, count = 0;
  for (k = 0; k < Txx.z; k++) {
    for (i = 0; i < Txx.x; i++) {
      fprintf(fp, "%le,", f[i][k]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}