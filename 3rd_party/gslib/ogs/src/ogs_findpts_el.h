#ifndef OGS_FINDPTS_EL_H
#define OGS_FINDPTS_EL_H

#if !defined(NAME_H) || !defined(POLY_H) || !defined(FINDPTS_EL_H)
#warning "ogs_findpts_el.h" requires "name.h", "poly.h", and "findpts_el.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void ogs_findpts_el_eval_2(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride, const unsigned pn,
  const double *const in, struct findpts_el_data_2 *const fd);

void ogs_findpts_el_eval_3(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride, const unsigned pn,
  const double *const in, struct findpts_el_data_3 *const fd);

#ifdef __cplusplus
}
#endif

#endif
