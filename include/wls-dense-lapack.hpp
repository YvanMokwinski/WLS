#pragma once
#pragma once
#include "wls-types.h"
#include "wls-mkl.hpp"

namespace WLS
{
  namespace dense
  {
    //!
    //! @brief BLAS MKL wrapper.
    //!
  
    struct lapack_t
    { 
      using int_t 	= wls_int_t;
      using int_p 	= wls_int_p;
      using const_int_p 	= const_wls_int_p;

      template <typename T> inline static void gesv(const_int_p 	m_,
						    const_int_p 	n_,
						    T*		x_,
						    const_int_p 	ld_,
						    wls_int_p 	ipiv_,
						    T*		y_,
						    const_int_p 	incy_,						     
						    wls_int_p	info_lapack_);

      template <typename T> inline static void getrs(const char * 	transpose_,
						     const_int_p 	m_,
						     const_int_p 	n_,
						     const T*	x_,
						     const_int_p 	ld_,
						     const_int_p 	ipiv_,
						     T*		y_,
						     const_int_p 	incy_,
						     wls_int_p	info_lapack_);

      template <typename T> inline static void getrf(const_int_p 	m_,
						     const_int_p 	n_,
						     T*	x_,
						     const_int_p 	ld_,
						     wls_int_p 	ipiv_,
						     wls_int_p		info_lapack_);


    };
  };
  

};
  
