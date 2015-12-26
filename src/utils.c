
#include "stdio.h"
#include "stdlib.h"
#include "utils.h"

int
make_linear(double * arr, int N, double a, double b)
{
  double h = (b - a) / (N - 1.);

  int I;
  FOREACH(I, N) arr[I] = a + h * I;

  return 0;
}

int
make_grid(double * ax_i, double * ax_j,
    int N_i, int N_j, double * X, double * Y)
{
  int I_i, I_j;

  for(I_i = 0; I_i < N_i; ++ I_i)
    for(I_j = 0; I_j < N_j; ++ I_j)
      X[I_i * N_j + I_j] = ax_i[I_i],
      Y[I_i * N_j + I_j] = ax_j[I_j];

  return 0;
}

int
save_array(double * arr, int size, const char * name)
{
  FILE * fp = fopen(name, "wb");
  if(fp == NULL) return -1;
  fwrite(arr, sizeof(double), size, fp);
  fclose(fp);
  return 0;
}

int
print_head(double * arr, int n)
{
  int I;
  FOREACH(I, n) printf("%f  ", arr[I]);
  printf("\n");
}

int
print_arr(double * arr, int n)
{
  int I;
  FOREACH(I, n) printf("%u %e\n", I, arr[I]);
  printf("\n");
}

int
draw_image(const char * name, int N_i, int N_j)
{
  char cmd_line[128];
  sprintf(cmd_line, "python debug-plot-2dimg.py %s %u %u", name, N_i, N_j);
  system(cmd_line);
  return 0;
}

int
draw_image_log(const char * name, int N_i, int N_j)
{
  char cmd_line[128];
  sprintf(cmd_line, "python debug-plot-2dimg-log.py %s %u %u", name, N_i, N_j);
  system(cmd_line);
  return 0;
}
