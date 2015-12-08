
#ifndef _SYNTHESIS_H
#define _SYNTHESIS_H

spectrum * synthesis(model * m_t, spec_lib * lib_t, double x, double y);
int synthesis_noalloc(model * m_t, spec_lib * lib_t,
    double x, double y, spectrum * sp_t);

#endif // _SYNTHESIS_H
