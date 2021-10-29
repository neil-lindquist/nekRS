/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

/* compile with C compiler (not C++) */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"
#include "ogs_findpts.h"

#include "ogstypes.h"

// need to access internals of findpts_data structs
struct hash_data_2 {
  ulong hash_n;
  struct dbl_range bnd[2];
  double fac[2];
  uint *offset;
};
struct findpts_data_2 {
  struct crystal cr;
  struct findpts_local_data_2 local;
  struct hash_data_2 hash;
};
struct hash_data_3 {
  ulong hash_n;
  struct dbl_range bnd[3];
  double fac[3];
  uint *offset;
};
struct findpts_data_3 {
  struct crystal cr;
  struct findpts_local_data_3 local;
  struct hash_data_3 hash;
};

void ogsDevFindpts_2(      dlong  *const  code_base   , const dlong  code_stride   ,
                           dlong  *const  proc_base   , const dlong  proc_stride   ,
                           dlong  *const    el_base   , const dlong    el_stride   ,
                           dfloat *const     r_base   , const dlong     r_stride   ,
                           dfloat *const dist2_base   , const dlong dist2_stride   ,
                     const dfloat *const     x_base[2], const dlong     x_stride[2],
                     const dlong npt, struct findpts_data_2 *const fd,
                     const void *const ogs_fd) {
  ogs_findpts_2( code_base,  code_stride,
                 proc_base,  proc_stride,
                   el_base,    el_stride,
                    r_base,     r_stride,
                dist2_base, dist2_stride,
                    x_base,     x_stride,
                npt, fd, ogs_fd);
}

void ogsDevFindpts_3(      dlong  *const  code_base   , const dlong  code_stride   ,
                           dlong  *const  proc_base   , const dlong  proc_stride   ,
                           dlong  *const    el_base   , const dlong    el_stride   ,
                           dfloat *const     r_base   , const dlong     r_stride   ,
                           dfloat *const dist2_base   , const dlong dist2_stride   ,
                     const dfloat *const     x_base[3], const dlong     x_stride[3],
                     const dlong npt, struct findpts_data_3 *const fd,
                     const void *const ogs_fd) {
  ogs_findpts_3( code_base,  code_stride,
                 proc_base,  proc_stride,
                   el_base,    el_stride,
                    r_base,     r_stride,
                dist2_base, dist2_stride,
                    x_base,     x_stride,
                npt, fd, ogs_fd);
}

void ogsDevFindptsEval_2(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_2 *const fd,
  const void *const ogs_fd) {

  ogs_findpts_eval_2( out_base,  out_stride,
                     code_base, code_stride,
                     proc_base, proc_stride,
                       el_base,   el_stride,
                        r_base,    r_stride,
                     npt, in, fd, ogs_fd);
}

void ogsDevFindptsEval_3(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_3 *const fd,
  const void *const ogs_fd) {

  ogs_findpts_eval_3( out_base,  out_stride,
                     code_base, code_stride,
                     proc_base, proc_stride,
                       el_base,   el_stride,
                        r_base,    r_stride,
                     npt, in, fd, ogs_fd);
}

void ogsDevFindptsLocalEval_2(
        void *const  out_base, const dlong  out_stride,
  const void *const   el_base, const dlong   el_stride,
  const void *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_2 *const fd,
  const void *const ogs_fd) {

  ogs_findpts_local_eval_2(out_base,  out_stride,
                            el_base,   el_stride,
                             r_base,    r_stride,
                           npt, in, &fd->local, ogs_fd);
}

void ogsDevFindptsLocalEval_3(
        void *const  out_base, const dlong  out_stride,
  const void *const   el_base, const dlong   el_stride,
  const void *const    r_base, const dlong    r_stride,
  const dlong npt, void *const in, struct findpts_data_3 *const fd,
  const void *const ogs_fd) {

  ogs_findpts_local_eval_3(out_base,  out_stride,
                            el_base,   el_stride,
                             r_base,    r_stride,
                           npt, in, &fd->local, ogs_fd);
}
