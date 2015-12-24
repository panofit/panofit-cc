
#include "../src/cps/sersic.h"
#include "../src/component.h"
#include "../src/model.h"
#include "../src/utils.h"
#include "../src/spectrum.h"
#include "../src/recipe.h"

#include "stdio.h"
#include "stdlib.h"

/*
#define IRS (par[ 0]) // surface density at the effective radius
#define XC  (par[ 1]) // center: x
#define YC  (par[ 2]) // center: y
#define PHI (par[ 3]) // pitch angle
#define RS  (par[ 4]) // effective radius
#define SN  (par[ 5]) // sersic n
#define Q   (par[ 6]) // axis ratio b/a
#define C   (par[ 7]) // shape parameter C
#define AMC (par[ 8])
#define ASC (par[ 9])
#define ZMC (par[10])
#define ZSC (par[11])
#define AMK (par[12])
#define ASK (par[13])
#define ZMK (par[14])
#define ZSK (par[15])
*/

int
main(int argc, char ** argv)
{
  // make and free sersic component // passed
  if(0)
    {
      double cp_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                         8., 0.2, 0., 0.2, 0., 0., 0., 0.};
      double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                         0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                         0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
      int cp_var[]    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

      component * cp_t =
        cps_sersic_2d(cp_par, cp_lim, cp_var, "bulge");

      free_component(cp_t);
    }

  // draw 2d surface brightness // passed
  if(0)
  {
    int I_pt;

    double cp_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                       8., 0.2, 0., 0.2, 0., 0., 0., 0.};
    double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                       0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                       0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    int cp_var[]    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    component * cp_t =
      cps_sersic_2d(cp_par, cp_lim, cp_var, "bulge");

    #define N_GRID_PTS 257

    double * ax = TALLOC(double, N_GRID_PTS),
           * Xt = TALLOC(double, N_GRID_PTS * N_GRID_PTS),
           * Yt = TALLOC(double, N_GRID_PTS * N_GRID_PTS),
           * It = TALLOC(double, N_GRID_PTS * N_GRID_PTS);

    make_linear(ax, N_GRID_PTS, -8., 8.);
    make_grid(ax, ax, N_GRID_PTS, N_GRID_PTS, Xt, Yt);

    _cps_sersic_2d_sigma_arr(N_GRID_PTS * N_GRID_PTS, Xt, Yt, cp_par, It);

    save_array(It, N_GRID_PTS * N_GRID_PTS, "sersic.tmp");
    draw_image_log("sersic.tmp", N_GRID_PTS, N_GRID_PTS);

    free_component(cp_t);
    free(ax), free(Xt), free(Yt), free(It);
  }

  // make recipe
  if(1)
  {
    int I_pt;

    double cp_par[] = {1., 0., 0., 0., 1., 2., 0.5, 2.,
                       8., 0.2, 0., 0.2, 0., 0., 0., 0.};
    double cp_lim[] = {1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 2., 2.,
                       0.5, 0.5, 2., 2., 8., 8., 0.2, 0.2, 0., 0.,
                       0.2, 0.2, 0., 0., 0., 0., 0., 0., 0., 0.};
    int cp_var[]    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    component * cp_t =
      cps_sersic_2d(cp_par, cp_lim, cp_var, "Terrible sersic");
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
