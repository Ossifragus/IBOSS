//' Obtain the index of the sample extremes
//' 
//' @param rr the number of minumum and maximum values to take
//' @param rz the vector to take extreme values from
//' @param rdel the elements to exclude from rz


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <algorithm> 
#include <iostream>

extern "C" {
  SEXP getIdxR(SEXP rr, SEXP rz, SEXP rdel) {
    SEXP idx;
    int r, *del, n, m;
    double *z;
    r = INTEGER(rr)[0];
    z = REAL(rz);
    n = length(rz);
    del = INTEGER(rdel);
    m = length(rdel);
    std::sort( del, del + m);
    double* y = new double [n-m];
    int j = 0, k = 0;
    for ( int i = 0; i < n; i++) {
      if ( j >= m) {
    	y[k++] = z[i];
      }
      else if ( del[j] != i + 1) {
    	y[k++] = z[i];
      }
      else
    	j++;
    }
    std::nth_element(y, y + r - 1, y + n - m);
    double  yrl = y[r-1];
    j = 0, k = 0;
    for ( int i = 0; i < n; i++) {
      if ( j >= m) {
    	y[k++] = -z[i];
      }
      else if ( del[j] != i + 1) {
    	y[k++] = -z[i];
      }
      else
    	j++;
    }
    std::nth_element(y, y + r - 1, y + n - m);
    double yru = -y[r-1];
    delete [] y;
    int jl = 0, ju = 0, locl[r], locu[r];
    j = 0;
    for ( int i = 0; i < n; i++) {
      if ( j >= m) {
	if ( z[i] <= yrl && jl < r)
	  locl[jl++] = i + 1;
	if ( z[i] >= yru && ju < r)
	  locu[ju++] = i + 1;
      }
      else if ( del[j] != i + 1) {
	if ( z[i] <= yrl && jl < r)
	  locl[jl++] = i + 1;
	if ( z[i] >= yru && ju < r)
	  locu[ju++] = i + 1;
      }
      else 
	j++;
      if ( jl >= r && ju >= r)
	break;
    }
    PROTECT(idx = allocVector(INTSXP, 2 * r));
    for ( int i = 0; i < r; i++) {
      INTEGER(idx)[i] = locl[i];
      INTEGER(idx)[r+i] = locu[i];
    }
    UNPROTECT(1);
    return idx;
  }


  SEXP getIdx(SEXP rr, SEXP rz) {
    SEXP idx;
    int r, n;
    double *z;
    r = INTEGER(rr)[0];
    z = REAL(rz);
    n = length(rz);
    double* y = new double [n];
    for ( int i = 0; i < n; i++)
      y[i] = z[i];
    std::nth_element(y, y + r - 1, y + n);
    double yrl = y[r-1];
    for ( int i = 0; i < n; i++)
      y[i] = -z[i];
    std::nth_element(y, y + r - 1, y + n);
    double yru = -y[r-1];
    delete [] y;
    int jl = 0, ju = 0, locl[r], locu[r];
    for ( int i = 0; i < n; i++) {
      if ( z[i] <= yrl && jl < r)
	locl[jl++] = i + 1;
      if ( z[i] >= yru && ju < r)
	locu[ju++] = i + 1;
      if ( jl >= r && ju >= r)
	break;
    }
    PROTECT(idx = allocVector(INTSXP, 2 * r));
    for ( int i = 0; i < r; i++) {
      INTEGER(idx)[i] = locl[i];
      INTEGER(idx)[r+i] = locu[i];
    }
    UNPROTECT(1);
    return idx;
  }
}

