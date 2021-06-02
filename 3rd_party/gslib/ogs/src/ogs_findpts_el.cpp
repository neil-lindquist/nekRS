
// work around the lack of restrict in C++
#define restrict

#include <cassert>
#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
extern "C" {
#include "name.h"
#include "poly.h"
#include "findpts_el.h"
}
#include "platform.hpp"
#include "ogs_findpts_el.h"

extern "C" {

void ogs_findpts_el_eval_2(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride,
  const unsigned pn, const double *const in, struct findpts_el_data_2 *const fd)
{
  const unsigned nr=fd->n[0],ns=fd->n[1];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);
  assert(nt <= MAX_GLL_N);

  assert(out_stride % sizeof(double) == 0);
  const unsigned d_out_stride = out_stride / sizeof(double);
  assert(r_stride % sizeof(double) == 0);
  const unsigned d_r_stride = r_stride / sizeof(double);

  occa::memory d_out_base = device.malloc(out_stride*pn);
  d_out_base.copyFrom(out_base);
  occa::memory   d_r_base = device.malloc(  r_stride*pn);
  d_r_base.copyFrom(r_base);
  occa::memory d_in = device.malloc(nr*ns*sizeof(double));
  d_in.copyFrom(in);

  occa::memory d_lag_data_0 = device.malloc(gll_lag_size(nr)*sizeof(double));
  d_lag_data_0.copyFrom(fd->lag_data[0]);
  occa::memory d_lag_data_1 = device.malloc(gll_lag_size(ns)*sizeof(double));
  d_lag_data_1.copyFrom(fd->lag_data[1]);

  ogs::findpts_el_eval_2(d_out_base, d_out_stride,
                           d_r_base,   d_r_stride,
                         pn, d_in,
                         nr, ns, d_lag_data_0, d_lag_data_1);

  d_out_base.copyTo(out_base);
}



void ogs_findpts_el_eval_3(
        double *const out_base, const unsigned out_stride,
  const double *const   r_base, const unsigned   r_stride,
  const unsigned pn, const double *const in, struct findpts_el_data_3 *const fd)
{
  const unsigned nr=fd->n[0],ns=fd->n[1],nt=fd->n[2];
  device_t device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);
  assert(nt <= MAX_GLL_N);

  assert(out_stride % sizeof(double) == 0);
  const unsigned d_out_stride = out_stride / sizeof(double);
  assert(r_stride % sizeof(double) == 0);
  const unsigned d_r_stride = r_stride / sizeof(double);

  occa::memory d_out_base = device.malloc(out_stride*pn);
  d_out_base.copyFrom(out_base);
  occa::memory   d_r_base = device.malloc(  r_stride*pn);
  d_r_base.copyFrom(r_base);
  occa::memory d_in = device.malloc(nr*ns*nt*sizeof(double));
  d_in.copyFrom(in);

  occa::memory d_lag_data_0 = device.malloc(gll_lag_size(nr)*sizeof(double));
  d_lag_data_0.copyFrom(fd->lag_data[0]);
  occa::memory d_lag_data_1 = device.malloc(gll_lag_size(ns)*sizeof(double));
  d_lag_data_1.copyFrom(fd->lag_data[1]);
  occa::memory d_lag_data_2 = device.malloc(gll_lag_size(nt)*sizeof(double));
  d_lag_data_2.copyFrom(fd->lag_data[2]);

  ogs::findpts_el_eval_3(d_out_base, d_out_stride,
                           d_r_base,   d_r_stride,
                         pn, d_in,
                         nr, ns, nt, d_lag_data_0, d_lag_data_1, d_lag_data_2);

  d_out_base.copyTo(out_base);
}

}
