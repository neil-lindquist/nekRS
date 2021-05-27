
#define interp_setup       TOKEN_PASTE(interp_setup_,D)
#define interp_free        TOKEN_PASTE(interp_free_ ,D)
#define interp_nfld        TOKEN_PASTE(interp_nfld  ,D)

#define findpts_setup       TOKEN_PASTE(findpts_setup_,D)
#define findpts_free        TOKEN_PASTE(findpts_free_ ,D)
#define findpts             TOKEN_PASTE(findpts_      ,D)
#define findpts_eval        TOKEN_PASTE(findpts_eval_ ,D)

// input:
//   nrs   ... nekRS configuration data
//   tol   ... tolerance newton solve (use 0 for default)
//
// return:
//   pointer to interpolation handles
struct findpts_data* interp_setup(nrs_t *nrs, double tol, unsigned nelm) {

  if (tol < 5e-13) {
    tolin = 5e-13;
  }
  int npt_max = 128;
  int bb_tol = 0.01;

  unsigned nmsh = nrs->mesh->N;
  unsigned nelm = nrs->mesh->Nelements;

  struct comm gs_comm;
  comm_init(&gs_comm, rns->mesh->comm);

  // element geometry
  double * elx[D];
  elx[0] = nrs->mesh->x;
  elx[1] = nrs->mesh->y;
  WHEN_3D(elx[2] = nrs->mesh->z);

  // element dimensions
  unsigned n1[D];
  n1[0] = nrs->mesh->N;
  n1[1] = nrs->mesh->N;
  WHEN_3D(n1[2] = nrs->mesh->N);

  unsigned m1[D];
  m1[0] = 2*n1[0];
  m1[1] = 2*n1[1];
  WHEN_3D(m1[2] = 2*n1[2]);

  // used for # of cells in hash tables
  int hash_size = nelm*n1[0]*n1[1] WHEN_3D(*n1[2]);

  return findpts_setup(gs_comm, elx, n1, nelm, m1, bb_tol,
                       hash_size, hash_size, npt_max, tol);
}


void interp_free(struct findpts_data *handle) {
  findpts_free(handle);
}

// input:
//   fld            ... source field(s)
//   nfld           ... number of fields
//   xp,yp,zp       ... interpolation points dim(n,D)
//   n              ... number of points
//   iwk            ... integer working array to hold point location information - dim[3,nmax]
//   rwk            ... real working array to hold the point local information - dim[D+1,nmax]
//   nmax           ... leading dimension of iwk and rwk
//   if_located_pts ... wheather to locate interpolation points (proc,el,r,s,t)
//   handle         ... handle
//
// output:
//   out            ... interpolation value(s) dim [nfld,n]
void interp_nfld(double *fld, unsigned nfld,
                 double *x[D],
                 int n, int *iwk, double *rwk,
                 int nmax, bool if_locate_pts,
                 struct findpts_data *handle,
                 double *out) {

  assert(n <= nmax);

  int    *code  = iwk;
  int    *proc  = iwk+2*nmax;
  int    *el    = iwk+nmax;
  double *r     = rwk+nmax;
  double *dist2 = rwk;

  bool if_trans_out = false;

  unsigned nfail = 0;
  if (if_locate_pts) {
    findpts(code,  1,
            proc,  1,
            el,    1,
            r,     D,
            dist2, 1,
            x,     1,
            n, handle);

    for (int in = 0; in < n; ++in) {
      if (code[in] == 1) {
        if (dist2[in] > 10*tol) {
          nfail += 1;
          //if (nfail < 5) write(6,'(a,1p4e15.7)')     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',     xp(in),yp(in),zp(in),rwk(in,1)
        }
      } else if (code[in] == 2) {
        nfail += 1
        //if (nfail < 5) write(6,'(a,1p3e15.7)')        ' WARNING: point not within mesh xy[z]: !',        xp(in),yp(in),zp(in)
      }
    }
  }

  for (int ifld = 0; ifld < nfld; ++ifld) {
     int in_offset  = ifld*nrs->fieldOffset;
     int out_offset = ifld*n;
     int out_stride = 1;
     if (if_trans_out) {
       out_offset = ifld;
       out_stride = nfld;
     }
     findpts_eval(out+out_offset, out_stride,
                  code,           1,
                  proc,           1,
                  el,             1,
                  r,              D,
                  n, fld+in_offset, handle);
  }

//  nn(1) = iglsum(n,1)
//  nn(2) = iglsum(nfail,1)
//  if(nio.eq.0) then
//    if(nn(2).gt.0 .or. loglevel.gt.2) write(6,1) nn(1),nn(2)
//1     format('   total number of points = ',i12,/,'   failed = '
// &         ,i12,/,' done :: intp_nfld')
//  endif
}


#undef interp_setup
#undef interp_free
#undef interp_nfld

#undef findpts_setup
#undef findpts_free
#undef findpts
#undef findpts_eval
