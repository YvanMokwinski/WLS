#pragma once
#include "wls-types.h"


namespace WLS
{
  namespace dense
  {
    //!
    //! @brief BLAS wrapper.
    //!  
    struct blas_t
    {    
      using int_t 	= wls_int_t;
      using int_p 	= wls_int_p;
      using const_int_p 	= const_wls_int_p;
    
      template <typename T>
      static inline void 	axpy(const_int_p 	n_,
				     const T* 		da_,
				     const T*		dx_,
				     const_int_p 	incx_,
				     T*			dy_,
				     const_int_p 	incy_);
    
      template <typename T>
      static inline void 	copy(const_int_p 	n_,
				     const T*		dx_,
				     const_int_p 	incx_,
				     T*			dy_,
				     const_int_p 	incy_);
    
      template <typename T>
      static inline void 	scal(const_int_p 	n_,
				     const T* 		da_,
				     T*			dx_,
				     const_int_p 	incx_);
    
      template <typename T>
      static inline T 	nrm2(const_int_p 	n_,
			     const T*		dx_,
			     const_int_p 	incx_);
    
      template <typename T>
      static inline void 	gemm(const char * 	transa_,
				     const char * 	transb_,
				     const_int_p 	m_,
				     const_int_p 	n_,
				     const_int_p 	k_,
				     const T* 		alpha_,
				     const T*		a_,
				     const_int_p 	lda_,
				     const T*		b_,
				     const_int_p 	ldb_,
				     const T* 		beta_,
				     T*			c_,
				     const_int_p 	ldc_);
    
      template <typename T>
      static inline void gemv(const char * 	trans_,
			      const_int_p 	m_,
			      const_int_p 	n_,
			      const T* 		alpha_,
			      const T*		a_,
			      const_int_p 	lda_,
			      const T*		x_,
			      const_int_p 	incx_,
			      const T* 		beta_,
			      T*			y_,
			      const_int_p 	incy_);
    };

  
  };
};
  
