
#ifndef _UTILS_H
#define _UTILS_H

// useful macro defs
#define TALLOC(_DATA_TYPE, _BLOCK_NUM) ((_DATA_TYPE *) malloc(sizeof(_DATA_TYPE) * (_BLOCK_NUM)))
#define FOREACH(IDX, NUM) for(IDX = 0; IDX < NUM; ++ IDX)
#define PI 3.1415926535897932

int make_linear(double *, int, double, double);
int make_grid(double *, double *, int, int, double *, double *);
int save_array(double *, int, const char *);

// debug functions
int print_arr(double *, int);
int print_head(double *, int);
int draw_image(const char *, int, int);
int draw_image_log(const char *, int, int);


#endif // _UTILS_H
