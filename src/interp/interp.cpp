
#include <mpi.h>
#include "nrs.hpp"
#include "platform.hpp"
#include <vector>

#include "interp.hpp"

// input:
//   nrs   ... nekRS configuration data
//   tol   ... tolerance newton solve (use 0 for default)
//
// return:
//   pointer to interpolation handles
struct interp_data* interp_setup(nrs_t *nrs, double tol) {

  if (tol < 5e-13) {
    tol = 5e-13;
  }
  int npt_max = 128;
  int bb_tol = 0.01;

  mesh_t *mesh = nrs->meshV;

  dlong nmsh = mesh->N;
  dlong nelm = mesh->Nelements;
  dlong D = mesh->dim;

  // element geometry
  dfloat *elx[3] = {mesh->x, mesh->y, mesh->z};

  // element dimensions
  dlong n1[3] = {mesh->N, mesh->N, mesh->N};

  dlong m1[3] = {2*n1[0], 2*n1[1], 2*n1[2]};

  // used for # of cells in hash tables
  dlong hash_size = nelm*n1[0]*n1[1];
  if (D == 3) hash_size *= n1[2];

  MPI_Comm comm = platform_t::getInstance()->comm.mpiComm;

  void *findpts_handle = ogsFindptsSetup(comm, D, elx, n1, nelm, m1, bb_tol,
                                         hash_size, hash_size, npt_max, tol);

  struct interp_data *handle = new interp_data();
  handle->nrs = nrs;
  handle->tol = tol;
  handle->findpts = findpts_handle;
}


void interp_free(struct interp_data *handle) {
  ogsFindptsFree((ogs_findpts_t*)handle->findpts);
  delete handle;
}

// input:
//   fld            ... source field(s)
//   nfld           ... number of fields
//   x              ... interpolation points dim[n,D]
//   n              ... number of points
//   iwk            ... integer working array to hold point location information - dim[3,nmax]
//   rwk            ... real working array to hold the point local information - dim[D+1,nmax]
//   nmax           ... leading dimension of iwk and rwk
//   if_located_pts ... wheather to locate interpolation points (proc,el,r,s,t)
//   handle         ... handle
//
// output:
//   out            ... interpolation value(s) dim [nfld,n]
void interp_nfld(dfloat *fld, dlong nfld,
                 dfloat *x[],
                 dlong n, dlong *iwk, dfloat *rwk,
                 dlong nmax, bool if_locate_pts,
                 struct interp_data *handle,
                 dfloat *out) {

  assert(n <= nmax);

  dlong  *code  = iwk;
  dlong  *proc  = iwk+2*nmax;
  dlong  *el    = iwk+nmax;
  dfloat *r     = rwk+nmax;
  dfloat *dist2 = rwk;

  bool if_trans_out = false;
  dlong D = handle->nrs->meshV->dim;
  dlong x_stride[3] = {1, 1, 1};

  unsigned nfail = 0;
  if (if_locate_pts) {
    ogsFindpts(code,  1,
               proc,  1,
               el,    1,
               r,     D,
               dist2, 1,
               x,     x_stride,
               n, (ogs_findpts_t*)handle->findpts);

    for (int in = 0; in < n; ++in) {
      if (code[in] == 1) {
        if (dist2[in] > 10*handle->tol) {
          nfail += 1;
          //if (nfail < 5) write(6,'(a,1p4e15.7)')     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',     xp(in),yp(in),zp(in),rwk(in,1)
        }
      } else if (code[in] == 2) {
        nfail += 1;
        //if (nfail < 5) write(6,'(a,1p3e15.7)')        ' WARNING: point not within mesh xy[z]: !',        xp(in),yp(in),zp(in)
      }
    }
  }

  for (int ifld = 0; ifld < nfld; ++ifld) {
     dlong in_offset  = ifld*handle->nrs->fieldOffset;
     dlong out_offset = ifld*n;
     dlong out_stride = 1;
     if (if_trans_out) {
       out_offset = ifld;
       out_stride = nfld;
     }
     ogsFindptsEval(out+out_offset, out_stride,
                    code,           1,
                    proc,           1,
                    el,             1,
                    r,              D,
                    n, fld+in_offset, (ogs_findpts_t*)handle->findpts);
  }

//  nn(1) = iglsum(n,1)
//  nn(2) = iglsum(nfail,1)
//  if(nio.eq.0) then
//    if(nn(2).gt.0 .or. loglevel.gt.2) write(6,1) nn(1),nn(2)
//1     format('   total number of points = ',i12,/,'   failed = '
// &         ,i12,/,' done :: intp_nfld')
//  endif
}
