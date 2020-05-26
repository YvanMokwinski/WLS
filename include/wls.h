#ifndef __wls_h__
#define __wls_h__

#include "wls-types.h"

extern "C"
{

  wls_status_t wls_matrix_market_dimension(FILE * 		in_,
					   wls_int_t * 	nrows_,
					   wls_int_t * 	ncols_,
					   wls_int_t * 	nnz_);
  
  wls_status_t wls_matrix_market_nnz	(FILE * 	in_,
					 wls_int_t 	direction_,
					 wls_int_t 	nrows_,
					 wls_int_t 	ncols_,
					 wls_int_t 	nnz_,
					 wls_int_p  	count_);
  
  wls_status_t wls_matrix_market_read_dcs(FILE * 	in_,
					  const char * 	idxsys_,
					  wls_int_t 	direction_,
					  wls_int_t  	nrows_,
					  wls_int_t  	ncols_,
					  wls_int_t  	nnz_,
					  wls_int_p  	begin_,
					  wls_int_p  	idx_,
					  double*__restrict__ 	values_);
  
  wls_status_t wls_matrix_market_read_scs	(FILE * 	in_,
						 const char * 	idxsys_,
						 wls_int_t 	direction_,
						 wls_int_t  	nrows_,
						 wls_int_t  	ncols_,
						 wls_int_t  	nnz_,
						 wls_int_p 	begin_,
						 wls_int_p 	idx_,
						 float*__restrict__ 	values_);

};

#endif
