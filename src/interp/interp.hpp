#if !defined(nekrs_interp_hpp_)
#define nekrs_interp_hpp_

#include "nrs.hpp"

// input:
//   nrs   ... nekRS configuration data
//   tol   ... tolerance newton solve (use 0 for default)
//
// return:
//   pointer to interpolation handles
struct findpts_data* interp_setup_2(nrs_t *nrs, double tol, unsigned nelm);
struct findpts_data* interp_setup_3(nrs_t *nrs, double tol, unsigned nelm);

void interp_free_2(struct findpts_data *handle);
void interp_free_3(struct findpts_data *handle);

// input:
//   fld            ... source field(s)
//   nfld           ... number of fields
//   xp,yp,zp       ... interpolation points dim(n,D)
//   n              ... number of points
//   iwk            ... integer working array to hold point location information - dim(nmax,3)
//   rwk            ... real working array to hold the canonical space coordinates and distances - dim(nmax,ldim+1)
//   nmax           ... leading dimension of iwk and rwk
//   if_located_pts ... wheather to locate interpolation points (proc,el,r,s,t)
//   handle         ... handle
//
// output:
//   out            ... interpolation value(s) dim (n,nfld)
void interp_nfld_2(double *fld, unsigned nfld,
                   double *x[2],
                   int n, int *iwk, double *rwk,
                   int nmax, bool if_locate_pts,
                   struct findpts_data *handle,
                   double *out);

void interp_nfld_3(double *fld, unsigned nfld,
                   double *x[3],
                   int n, int *iwk, double *rwk,
                   int nmax, bool if_locate_pts,
                   struct findpts_data *handle,
                   double *out);


#endif
