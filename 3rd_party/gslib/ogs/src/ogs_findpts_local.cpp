
// work around the lack of restrict in C++
#define restrict

#include <cassert>
#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "platform.hpp"

#define   AT(T,var,i)   \
        (T*)(      (char*)var##_base   +(i)*var##_stride   )
#define  CAT(T,var,i) \
  (const T*)((const char*)var##_base   +(i)*var##_stride   )
#define CATD(T,var,i,d) \
  (const T*)((const char*)var##_base[d]+(i)*var##_stride[d])

extern "C" {

void ogs_findpts_local_eval_internal_2(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const unsigned pn, const double *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[2], unsigned lag_data_size[2])
{
  if (pn == 0) return;

  const unsigned nr=n[0],ns=n[1];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);

  assert(out_stride % sizeof(double) == 0);
  const unsigned d_out_stride = out_stride / sizeof(double);
  assert(r_stride % sizeof(double) == 0);
  const unsigned d_r_stride = r_stride / sizeof(double);
  assert(el_stride % sizeof(double) == 0);
  const unsigned d_el_stride = el_stride / sizeof(uint);

  unsigned max_el = 0;
  for (unsigned i = 0; i < pn; ++i) {
    if (max_el < *CAT(unsigned, el, i)) max_el = *CAT(unsigned, el, i);
  }

  occa::json malloc_props;

  occa::memory d_out_base = device.malloc(out_stride*pn,
                                          occa::dtype::byte,
                                          malloc_props);
  d_out_base.copyFrom(out_base);
  occa::memory d_el_base = device.malloc(el_stride*pn,
                                         occa::dtype::byte,
                                         malloc_props);
  d_el_base.copyFrom(el_base);
  occa::memory   d_r_base = device.malloc(r_stride*pn,
                                          occa::dtype::byte,
                                          malloc_props);
  d_r_base.copyFrom(r_base);

  occa::memory d_in = device.malloc(in_stride*(max_el+1), occa::dtype::double_, malloc_props);
  d_in.copyFrom(in);

  occa::memory d_lag_data_0 = device.malloc(lag_data_size[0],
                                            occa::dtype::double_,
                                            malloc_props);
  d_lag_data_0.copyFrom(lag_data[0]);
  // reuse lag_data_0 if all directions have the same degree
  occa::memory d_lag_data_1;
  if (nr == ns) {
    d_lag_data_1 = d_lag_data_0;
  } else {
    d_lag_data_1 = device.malloc(lag_data_size[1],
                                            occa::dtype::double_,
                                            malloc_props);
    d_lag_data_1.copyFrom(lag_data[1]);
  }

  ogs::findpts_local_eval_2(d_out_base, d_out_stride,
                             d_el_base,  d_el_stride,
                              d_r_base,   d_r_stride,
                            pn, d_in, in_stride,
                            nr, ns, d_lag_data_0, d_lag_data_1);

  d_out_base.copyTo(out_base);
}



void ogs_findpts_local_eval_internal_3(
        double *const out_base, const unsigned out_stride,
  const uint   *const  el_base, const unsigned  el_stride,
  const double *const   r_base, const unsigned   r_stride,
  const unsigned pn, const double *const in, const unsigned in_stride,
  unsigned *const n, double *const lag_data[3], unsigned lag_data_size[3])
{
  if (pn == 0) return;

  const unsigned nr=n[0],ns=n[1],nt=n[2];
  occa::device device = platform->device;

  assert(nr <= MAX_GLL_N);
  assert(ns <= MAX_GLL_N);
  assert(nt <= MAX_GLL_N);

  assert(out_stride % sizeof(double) == 0);
  const unsigned d_out_stride = out_stride / sizeof(double);
  assert(r_stride % sizeof(double) == 0);
  const unsigned d_r_stride = r_stride / sizeof(double);
  assert(el_stride % sizeof(double) == 0);
  const unsigned d_el_stride = el_stride / sizeof(uint);

  unsigned max_el = 0;
  for (unsigned i = 0; i < pn; ++i) {
    if (max_el < *CAT(unsigned, el, i)) max_el = *CAT(unsigned, el, i);
  }

  occa::json malloc_props;

  occa::memory d_out_base = device.malloc(out_stride*pn,
                                          occa::dtype::byte,
                                          malloc_props);
  d_out_base.copyFrom(out_base);
  occa::memory d_el_base = device.malloc(el_stride*pn,
                                         occa::dtype::byte,
                                         malloc_props);
  d_el_base.copyFrom(el_base);
  occa::memory   d_r_base = device.malloc(r_stride*pn,
                                          occa::dtype::byte,
                                          malloc_props);
  d_r_base.copyFrom(r_base);

  occa::memory d_in = device.malloc(in_stride*(max_el+1), occa::dtype::double_, malloc_props);
  d_in.copyFrom(in);

  occa::memory d_lag_data_0 = device.malloc(lag_data_size[0],
                                            occa::dtype::double_,
                                            malloc_props);
  d_lag_data_0.copyFrom(lag_data[0]);
  // reuse lag_data_0 if all directions have the same degree
  occa::memory d_lag_data_1, d_lag_data_2;
  if (nr == ns) {
    d_lag_data_1 = d_lag_data_0;
  } else {
    d_lag_data_1 = device.malloc(lag_data_size[1], occa::dtype::double_, malloc_props);
    d_lag_data_1.copyFrom(lag_data[1]);
  }
  if (nr == nt) {
    d_lag_data_2 = d_lag_data_0;
  } else {
    d_lag_data_2 = device.malloc(lag_data_size[2], occa::dtype::double_, malloc_props);
    d_lag_data_2.copyFrom(lag_data[2]);
  }

  ogs::findpts_local_eval_3(d_out_base, d_out_stride,
                             d_el_base,  d_el_stride,
                              d_r_base,   d_r_stride,
                            pn, d_in, in_stride,
                            nr, ns, nt, d_lag_data_0, d_lag_data_1, d_lag_data_2);

  d_out_base.copyTo(out_base);
}

}
