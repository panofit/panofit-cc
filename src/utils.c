
#include "model.h"

#include "stdio.h"
#include "stdlib.h"
#include "utils.h"

#include "math.h"
#include "float.h"

#define Sq(X) ((X) * (X))

double
rand_gauss()
{
  double        u, v, w, m;
  static double x, y;
  static int    p = 0;

  if (p == 1) {p = !p; return (double) y;}

  do
    {
      u = -1. + ((double) rand () / RAND_MAX) * 2.;
      v = -1. + ((double) rand () / RAND_MAX) * 2.;
      w = Sq(u) + Sq(v);
    }
  while (w >= 1. || w == 0.);

  m = sqrt ((-2. * log (w)) / w);
  x = u * m, y = v * m;
  p = !p;

  return (double) x;
}

double
rand_uniform()
{
  return ((double) rand () / RAND_MAX);
}

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
print_param_set(param_set * ps_t, int variable_only)
{
  // iterate over components
  int N_cps = ps_t -> N_cps, I_cp;
  FOREACH(I_cp, N_cps)
    {
      // print component name
      printf("Component %u/%u: %s\n", I_cp + 1, N_cps, ps_t -> name[I_cp]);

      // iterate over parameters, print names and values
      int I_par, N_par = ps_t -> N_par[I_cp];
      FOREACH(I_par, N_par)
        {
          if(!(I_par % 4)) printf("\n  "); // just formatting

          // print with yellow text if variable
          if(!((ps_t -> is_const[I_cp])[I_par])) // variable
            printf("\x1B[32m%s: %f\033[0m,  ",
                (ps_t -> par_name[I_cp])[I_par], (ps_t -> par[I_cp])[I_par]);
          else if(! variable_only)
            printf("%s: %f,  ",
                (ps_t -> par_name[I_cp])[I_par], (ps_t -> par[I_cp])[I_par]);
        }

      printf("\n\n");
    }
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
