
#include "../src/recipe.h"
#include "../src/spectrum.h"

#include "stdio.h"

int
main(int argc, char ** argv)
{
  // read spectral library and free // passed
  if(0)
    {
      int I_pt;

      // read spec lib
      spec_lib * lib_t = load_spec_lib_raw("mock-speclib.dat");

      // print name
      printf("Name: %s\n", lib_t -> name);
      printf("Age pts: %u, Meta. pts: %u, wl pts: %u\n",
          lib_t -> N_age, lib_t -> N_Z, lib_t -> N_spx);

      // print age axis
      printf("age axis: ");
      for(I_pt = 0; I_pt < lib_t -> N_age; ++ I_pt)
        printf("%.2f ", lib_t -> age_ax[I_pt]);
      printf("\n");

      // print metallicity axis
      printf("metallicity axis: ");
      for(I_pt = 0; I_pt < lib_t -> N_Z; ++ I_pt)
        printf("%.2f ", lib_t -> Z_ax[I_pt]);
      printf("\n");

      // print wavelength axis
      printf("wavelength axis: ");
      for(I_pt = 0; I_pt < lib_t -> N_spx; ++ I_pt)
        printf("%.2f ", lib_t -> wl[I_pt]);
      printf("\n");

      // print the first ssp
      printf("first ssp:\n");
      for(I_pt = 0; I_pt < lib_t -> N_spx; ++ I_pt)
        printf("%.2f, %e\n", lib_t -> wl[I_pt], lib_t -> data[I_pt]);
      printf("\n");

      free_spec_lib(lib_t);
    }

  // create empty spectrum // passed
  if(0)
    {
      // read spec lib
      spec_lib * lib_t = load_spec_lib_raw("mock-speclib.dat");

      // create empty spectrum
      spectrum * sp_t = make_empty_spectrum(lib_t);

      // print wavelength axis in sp_t
      printf("wavelength axis: ");
      int I_pt;
      for(I_pt = 0; I_pt < lib_t -> N_spx; ++ I_pt)
        printf("%.2f ", sp_t -> wl[I_pt]);
      printf("\n");

      // free spec.
      free_spectrum(sp_t);

      // print wavelength axis in lib_t
      printf("wavelength axis: ");
      for(I_pt = 0; I_pt < lib_t -> N_spx; ++ I_pt)
        printf("%.2f ", lib_t -> wl[I_pt]);
      printf("\n");

      // free spec_lib
      free_spec_lib(lib_t);
    }

    // create empty recipe // passed
    if(0)
      {
        int I_pt;

        // read spec lib
        spec_lib * lib_t = load_spec_lib_raw("mock-speclib.dat");

        // create empty spectrum
        recipe * rc_t = make_empty_recipe(lib_t);

        // print age axis
        printf("age axis: ");
        for(I_pt = 0; I_pt < rc_t -> N_age; ++ I_pt)
          printf("%.2f ", rc_t -> age_ax[I_pt]);
        printf("\n");

        // print metallicity axis
        printf("metallicity axis: ");
        for(I_pt = 0; I_pt < rc_t -> N_Z; ++ I_pt)
          printf("%.2f ", rc_t -> Z_ax[I_pt]);
        printf("\n");

        // free recipe and spec lib
        free_recipe(rc_t), free_spec_lib(lib_t);
      }

}
