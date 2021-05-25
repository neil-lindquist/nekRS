#ifndef FINDPTS_EL_H
#define FINDPTS_EL_H

#if !defined(NAME_H) || !defined(POLY_H)
#warning "ogs_findpts_el.h" requires "poly.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


struct findpts_el_pt_2 {
  double x[2],r[2],oldr[2],dist2,dist2p,tr;
  unsigned index,flags;
};

struct findpts_el_gedge_2 { const double *x[2], *dxdn[2]; };
struct findpts_el_gpt_2   { double x[2], jac[4], hes[4]; };

struct findpts_el_data_2 {
  unsigned npt_max;
  struct findpts_el_pt_2 *p;

  unsigned n[2];
  double *z[2];
  lagrange_fun *lag[2];
  double *lag_data[2];
  double *wtend[2];

  const double *x[2];

  unsigned side_init;
  double *sides;
  struct findpts_el_gedge_2 edge[4]; /* R S=-1; R S=1; ... */
  struct findpts_el_gpt_2 pt[4];

  double *work;
};


struct findpts_el_pt_3 {
  double x[3],r[3],oldr[3],dist2,dist2p,tr;
  unsigned index,flags;
};

struct findpts_el_gface_3 { const double *x[3], *dxdn[3]; };
struct findpts_el_gedge_3 { const double *x[3], *dxdn1[3], *dxdn2[3],
                                         *d2xdn1[3], *d2xdn2[3]; };
struct findpts_el_gpt_3   { double x[3], jac[9], hes[18]; };

struct findpts_el_data_3 {
  unsigned npt_max;
  struct findpts_el_pt_3 *p;

  unsigned n[3];
  double *z[3];
  lagrange_fun *lag[3];
  double *lag_data[3];
  double *wtend[3];

  const double *x[3];

  unsigned side_init;
  double *sides;
  struct findpts_el_gface_3 face[6]; /* ST R=-1,R=+1; TR S=-1,S=+1; ... */
  struct findpts_el_gedge_3 edge[12]; /* R S=-1,T=-1; R S=1,T=-1; ... */
  struct findpts_el_gpt_3 pt[8];

  double *work;
};

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
