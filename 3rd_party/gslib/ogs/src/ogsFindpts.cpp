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

#include <cassert>
#include <cstdlib>
#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsInterface.h"

ogs_findpts_t *ogsFindptsSetup(
  const dlong D, MPI_Comm comm,
  const dfloat *const elx[],
  const dlong n[], const dlong nel,
  const dlong m[], const dfloat bbox_tol,
  const hlong local_hash_size, const hlong global_hash_size,
  const dlong npt_max, const dfloat newt_tol) {
  // elx, n, m have length D

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  void *findpts_data;
  if (D == 2) {
    findpts_data = ogsHostFindptsSetup_2(comm, elx, n, nel, m, bbox_tol,
                                         local_hash_size, global_hash_size,
                                         npt_max, newt_tol);
  } else if (D == 3) {
    findpts_data = ogsHostFindptsSetup_3(comm, elx, n, nel, m, bbox_tol,
                                         local_hash_size, global_hash_size,
                                         npt_max, newt_tol);
  } else {
    assert(false);
  }

  ogs_findpts_t *ogs_handle = (ogs_findpts_t*)malloc(sizeof(ogs_findpts_t));
  ogs_handle->D = D;
  ogs_handle->findpts_data = findpts_data;
  return ogs_handle;
}


void ogsFindptsFree(ogs_findpts_t *fd) {
  if (fd->D == 2) {
    ogsHostFindptsFree_2((findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindptsFree_3((findpts_data_3*)fd->findpts_data);
  }
  free(fd);
}

void ogsFindpts(    dlong  *const  code_base  , const dlong  code_stride,
                    dlong  *const  proc_base  , const dlong  proc_stride,
                    dlong  *const    el_base  , const dlong    el_stride,
                    dfloat *const     r_base  , const dlong     r_stride,
                    dfloat *const dist2_base  , const dlong dist2_stride,
              const dfloat *const     x_base[], const dlong     x_stride[],
              const dlong npt, ogs_findpts_t *const fd) {
  // x_base, x_stride have length D

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  if (fd->D == 2) {
    ogsHostFindpts_2( code_base,  code_stride,
                      proc_base,  proc_stride,
                        el_base,    el_stride,
                         r_base,     r_stride,
                     dist2_base, dist2_stride,
                         x_base,     x_stride,
                     npt, (findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindpts_3( code_base,  code_stride,
                      proc_base,  proc_stride,
                        el_base,    el_stride,
                         r_base,     r_stride,
                     dist2_base, dist2_stride,
                         x_base,     x_stride,
                     npt, (findpts_data_3*)fd->findpts_data);
  }
}

void ogsFindptsEval(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, const dfloat *const in, ogs_findpts_t *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  if (fd->D == 2) {
    ogsHostFindptsEval_2( out_base,  out_stride,
                         code_base, code_stride,
                         proc_base, proc_stride,
                           el_base,   el_stride,
                            r_base,    r_stride,
                         npt, in, (findpts_data_2*)fd->findpts_data);
  } else {
    ogsHostFindptsEval_3( out_base,  out_stride,
                         code_base, code_stride,
                         proc_base, proc_stride,
                           el_base,   el_stride,
                            r_base,    r_stride,
                         npt, in, (findpts_data_3*)fd->findpts_data);
  }
}


void ogsFindptsEval(
        dfloat *const  out_base, const dlong  out_stride,
  const dlong  *const code_base, const dlong code_stride,
  const dlong  *const proc_base, const dlong proc_stride,
  const dlong  *const   el_base, const dlong   el_stride,
  const dfloat *const    r_base, const dlong    r_stride,
  const dlong npt, occa::memory d_in, ogs_findpts_t *const fd) {

  assert(sizeof(dfloat) == sizeof(double));
  assert(sizeof(dlong) == sizeof(uint));

  if (fd->D == 2) {
    ogsDevFindptsEval_2( out_base,  out_stride,
                        code_base, code_stride,
                        proc_base, proc_stride,
                          el_base,   el_stride,
                           r_base,    r_stride,
                        npt, &d_in, (findpts_data_2*)fd->findpts_data);
  } else {
    ogsDevFindptsEval_3( out_base,  out_stride,
                        code_base, code_stride,
                        proc_base, proc_stride,
                          el_base,   el_stride,
                           r_base,    r_stride,
                        npt, &d_in, (findpts_data_3*)fd->findpts_data);
  }
}
