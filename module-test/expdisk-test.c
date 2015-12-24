
#include "../src/cps/expdisk.h"
#include "../src/component.h"
#include "../src/model.h"
#include "../src/utils.h"
#include "../src/spectrum.h"
#include "../src/recipe.h"

#include "stdio.h"
#include "stdlib.h"

/*
#define IC  (par[ 0])
#define XC  (par[ 1])
#define YC  (par[ 2])
#define PHI (par[ 3])
#define A   (par[ 4])
#define Q   (par[ 5])
#define C   (par[ 6])
#define AMC (par[ 7])
#define ASC (par[ 8])
#define ZMC (par[ 9])
#define ZSC (par[10])
#define AMK (par[11])
#define ASK (par[12])
#define ZMK (par[13])
#define ZSK (par[14])
*/

int
main(int argc, char ** argv)
{
  // make and free exponential disk // passed
  if(0)
    {
      double cp_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                         8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                         2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                         0., 0., 0., 0., 0., 0.};
      int cp_var[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

      component * cp_t =
        cps_expdisk_2d(cp_par, cp_lim, cp_var, "Lovely expdisk");

      free_component(cp_t);
    }

  // make 2d image // passed
  if(0)
  {
    int I_pt;

    double cp_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                       8., 0.2, 0., 0.2, 0., 0., 0., 0.};
    double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                       2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                       0., 0., 0., 0., 0., 0.};
    int cp_var[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    component * cp_t =
      cps_expdisk_2d(cp_par, cp_lim, cp_var, "Lovely expdisk");

    #define N_GRID_PTS 257

    double * ax = TALLOC(double, N_GRID_PTS),
           * Xt = TALLOC(double, N_GRID_PTS * N_GRID_PTS),
           * Yt = TALLOC(double, N_GRID_PTS * N_GRID_PTS),
           * It = TALLOC(double, N_GRID_PTS * N_GRID_PTS);

    make_linear(ax, N_GRID_PTS, -8., 8.);
    make_grid(ax, ax, N_GRID_PTS, N_GRID_PTS, Xt, Yt);

    _cps_expdisk_2d_sigma_arr(N_GRID_PTS * N_GRID_PTS, Xt, Yt, cp_par, It);

    save_array(It, N_GRID_PTS * N_GRID_PTS, "expdisk.tmp");
    draw_image_log("expdisk.tmp", N_GRID_PTS, N_GRID_PTS);

    free_component(cp_t);
    free(ax), free(Xt), free(Yt), free(It);
  }

  // make recipe // passed
  if(0)
  {
    int I_pt;

    double cp_par[] = {1., 0., 0., 0., 1., 0.5, 2.,
                       10., 1., 0., 0.2, 0., 0., 0., 0.};
    double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 0.5, 0.5, 2.,
                       2., 8., 8., 0.2, 0.2, 0., 0., 0.2, 0.2, 0., 0.,
                       0., 0., 0., 0., 0., 0.};
    int    cp_var[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    component * cp_t =
      cps_expdisk_2d(cp_par, cp_lim, cp_var, "Lovely expdisk");
    model * m_t = make_model("Playful model", 1, cp_t);

    // load spectral library
    spec_lib * lib_t = load_spec_lib_raw("mock-speclib.dat");

    // create empty recipe
    recipe * rc_t = make_empty_recipe(lib_t);

    // generate recipe
    sample_recipe_noalloc(m_t, 1., 1., rc_t);

    save_array(rc_t -> rcp, (rc_t -> N_Z) * (rc_t -> N_age), "rcp.tmp");
    draw_image("rcp.tmp", (rc_t -> N_age), (rc_t -> N_Z));
  }
}
