
#ifndef _UTILS_H
#define _UTILS_H

#define TALLOC(type, N) ((type *) malloc(sizeof(type) * (N)))
#define FOREACH(IDX, NUM) for(IDX = 0; IDX < NUM; ++ IDX)

#endif // _UTILS_H
