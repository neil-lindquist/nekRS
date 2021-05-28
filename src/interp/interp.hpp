#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nrs.hpp"

struct interp_data {
  nrs_t *nrs;
  double tol;
  unsigned D;
  void *findpts;
};

// input:
//   nrs   ... nekRS configuration data
//   tol   ... tolerance newton solve (use 0 for default)
//
// return:
//   pointer to interpolation handles
struct interp_data* interp_setup(nrs_t *nrs, double tol);
void interp_free(struct interp_data *handle);

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
                 dfloat *x[], dlong x_stride[],
                 dlong n, dlong *iwk, dfloat *rwk,
                 dlong nmax, bool if_locate_pts,
                 struct interp_data *handle,
                 dfloat *out, bool if_trans_out = false);


#endif
