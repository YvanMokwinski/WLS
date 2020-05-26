#pragma once

#include "wls-dense-lapack.hpp"

namespace WLS
{
  namespace dense
  {

    template <> inline void lapack_t::gesv<float>(const_int_p 	m_,
						  const_int_p 	n_,
						  float*		x_,
						  const_int_p 	ld_,
						  wls_int_p 	ipiv_,
						  float*		y_,
						  const_int_p 	incy_,						     
						  wls_int_p	info_lapack_)
    {
    
      sgesv(m_,
	    n_,
	    x_,
	    ld_,
	    ipiv_,
	    y_,
	    incy_,
	    info_lapack_);

    };

    template <> inline void lapack_t::gesv<double>(const_int_p 	m_,
						   const_int_p 	n_,
						   double*		x_,
						   const_int_p 	ld_,
						   wls_int_p 	ipiv_,
						   double*		y_,
						   const_int_p 	incy_,						     
						   wls_int_p	info_lapack_)
    {

      dgesv(m_,
	    n_,
	    x_,
	    ld_,
	    ipiv_,
	    y_,
	    incy_,
	    info_lapack_);

    };

    template <> inline void lapack_t::getrs<float>(const char * 	transpose_,
						   const_int_p 	m_,
						   const_int_p 	n_,
						   const float*	x_,
						   const_int_p 	ld_,
						   const_int_p 	ipiv_,
						   float*		y_,
						   const_int_p 	incy_,
						   wls_int_p	info_lapack_)
    {

      sgetrs(transpose_,
	     m_,
	     n_,
	     x_,
	     ld_,
	     ipiv_,
	     y_,
	     incy_,
	     info_lapack_);

    };

    template <> inline void lapack_t::getrs<double>(const char * 	transpose_,
						    const_int_p 	m_,
						    const_int_p 	n_,
						    const double*	x_,
						    const_int_p 	ld_,
						    const_int_p 	ipiv_,
						    double*		y_,
						    const_int_p 	incy_,
						    wls_int_p	info_lapack_)
    {

      dgetrs(transpose_,
	     m_,
	     n_,
	     x_,
	     ld_,
	     ipiv_,
	     y_,
	     incy_,
	     info_lapack_);

    };

    template <> inline void lapack_t::getrf<float>(const_int_p 	m_,
						   const_int_p 	n_,
						   float*	x_,
						   const_int_p 	ld_,
						   wls_int_p 	ipiv_,
						   wls_int_p		info_lapack_)
    {

      sgetrf(m_,
	     n_,
	     x_,
	     ld_,
	     ipiv_,
	     info_lapack_);

    };

    template <> inline void lapack_t::getrf<double>(const_int_p 	m_,
						    const_int_p 	n_,
						    double*	x_,
						    const_int_p 	ld_,
						    wls_int_p 	ipiv_,
						    wls_int_p		info_lapack_)
    {

      dgetrf(m_,
	     n_,
	     x_,
	     ld_,
	     ipiv_,
	     info_lapack_);

    };
  };
};
  
