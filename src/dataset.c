
#include "dataset.h"
#include "utils.h"

dataset *
make_dataset(const char * name, int N_img, int N_spec)
{
  // create empty dataset
  dataset * ds_t = TALLOC(dataset, 1);
  ds_t -> img = TALLOC(image *, N_img),
  ds_t -> spec = TALLOC(spec *, N_spec);

  // initialize (set null)
  for I_pc;
  FOREACH(I_pc, N_img) ds_t -> img[I_pc] = NULL;
  FOREACH(I_pc, N_spec) ds_t -> spec[I_pc] = NULL;

  // write properties
  ds_t -> N_img = N_img,
  ds_t -> N_spec = N_spec;
  strcpy(ds_t -> name, name);

  // return object
  return ds_t;
}

int
free_dataset(dataset * ds_t, int free_objs)
{
  if(free_objs)
    {
      int I_pc;

      // delete image objects
      FOREACH(I_pc, ds_t -> N_img)
        if(ds_t -> img[I_pc] != NULL) free(ds_t -> img[I_pc]);

      // delete spectra
      FOREACH(I_pc, ds_t -> N_spec)
        if(ds_t -> spec[I_pc] != NULL) free(ds_t -> spec[I_pc]);
    }

  free(ds_t -> img), free(ds_t -> spec);
  free(ds_t); ds_t = NULL;

  return 0;
}
