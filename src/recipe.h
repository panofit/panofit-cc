
#ifndef _RECIPE_H
#define _RECIPE_H

typedef struct _recipe
{
  double x, y, z; // coord in the observer's frame

  double * rcp;

} recipe;

recipe make_empty_recipe();

#endif // _RECIPE_H
