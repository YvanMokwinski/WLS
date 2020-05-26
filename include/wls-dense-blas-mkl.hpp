#pragma once
#include "wls-dense-blas.hpp"
#include "wls-mkl.hpp"

namespace WLS
{
  namespace dense
  {
    
    template <>
    inline float 	blas_t::nrm2<float>	(const_int_p 	n_,
							 const float*	dx_,
							 const_int_p 	incx_)
    {
    
      return snrm2(n_,
		   dx_,
		   incx_);
    
    };
  
    template <>
    inline double 	blas_t::nrm2<double>	(const_int_p 	n_,
							 const double*	dx_,
							 const_int_p 	incx_)
    {
    
      return dnrm2(n_,
		   dx_,
		   incx_);
    };
  
    template <>
    inline void 	blas_t::scal<float>	(const_int_p 	n_,
						 const float* 	da_,
						 float*		dx_,
						 const_int_p 	incx_)
    {
    
      sscal((wls_int_p)n_,
	    (float*)da_,
	    (float*)dx_,
	    (wls_int_p)incx_);
    };
  
    template <>
    inline void 	blas_t::scal<double>	(const_int_p 	n_,
						 const double* 	da_,
						 double*	dx_,
						 const_int_p 	incx_)
    {
    
      dscal((wls_int_p)n_,
	    (double*)da_,
	    (double*)dx_,
	    (wls_int_p)incx_);
    };
  
    template <> inline void blas_t::axpy<float>	(const_int_p 	n_,
						 const float* 	da_,
						 const float*	dx_,
						 const_int_p 	incx_,
						 float*		dy_,
						 const_int_p 	incy_)
    {
    
      saxpy((wls_int_p)n_,
	    (float*)da_,
	    (float*)dx_,
	    (wls_int_p)incx_,
	    (float*)dy_,
	    (wls_int_p)incy_);
    
    
    };
  
    template <> inline void blas_t::axpy<double>(const_int_p 	n_,
						 const double* 			da_,
						 const double*			dx_,
						 const_int_p 	incx_,
						 double*				dy_,
						 const_int_p 	incy_)
    {
    
      daxpy((wls_int_p)n_,
	    (double*)da_,
	    (double*)dx_,
	    (wls_int_p)incx_,
	    (double*)dy_,
	    (wls_int_p)incy_);
    
    };
  
    template <> inline void blas_t::copy<float>(const_int_p 	n_,
						const float*			dx_,
						const_int_p 	incx_,
						float*				dy_,
						const_int_p 	incy_)
    {
      scopy((wls_int_p)n_,
	    (float*)dx_,
	    (wls_int_p)incx_,
	    (float*)dy_,
	    (wls_int_p)incy_);

    };

    template <> inline void blas_t::copy<double>(const_int_p 	n_,
						 const double*			dx_,
						 const_int_p 	incx_,
						 double*				dy_,
						 const_int_p 	incy_)
    {

      dcopy((wls_int_p)n_,
	    (double*)dx_,
	    (wls_int_p)incx_,
	    (double*)dy_,
	    (wls_int_p)incy_);

    };


    template <> inline void blas_t::gemm<double>(const char * 			transa_,
						 const char * 			transb_,
						 const_int_p 	m_,
						 const_int_p 	n_,
						 const_int_p 	k_,
						 const double* 			alpha_,
						 const double*			a_,
						 const_int_p 	lda_,
						 const double*b_,
						 const_int_p ldb_,
						 const double* beta_,
						 double*c_,
						 const_int_p ldc_)
    {

      dgemm((char*)transa_,
	    (char*)transb_,
	    (wls_int_p)m_,
	    (wls_int_p)n_,
	    (wls_int_p)k_,
	    (double*)alpha_,
	    (double*)a_,
	    (wls_int_p)lda_,
	    (double*)b_,
	    (wls_int_p)ldb_,
	    (double*)beta_,
	    c_,
	    (wls_int_p)ldc_);

    };

    template <> inline void blas_t::gemm<float>(const char * 			transa_,
						const char * 			transb_,
						const_int_p 	m_,
						const_int_p 	n_,
						const_int_p 	k_,
						const float* 			alpha_,
						const float*			a_,
						const_int_p 	lda_,
						const float*			b_,
						const_int_p 	ldb_,
						const float* 			beta_,
						float*				c_,
						const_int_p 	ldc_)
    {

      sgemm((char*)transa_,
	    (char*)transb_,
	    (wls_int_p)m_,
	    (wls_int_p)n_,
	    (wls_int_p)k_,
	    (float*)alpha_,
	    (float*)a_,
	    (wls_int_p)lda_,
	    (float*)b_,
	    (wls_int_p)ldb_,
	    (float*)beta_,
	    c_,
	    (wls_int_p)ldc_);

    };
  

    template <> inline void blas_t::gemv<float>(const char * trans_,
						const_int_p m_,
						const_int_p n_,
						const float* alpha_,
						const float*a_,
						const_int_p lda_,
						const float*x_,
						const_int_p incx_,
						const float* beta_,
						float*y_,
						const_int_p incy_)
    {

      sgemv((char*)trans_,
	    (wls_int_p)m_,
	    (wls_int_p)n_,
	    (float*)alpha_,
	    (float*)a_,
	    (wls_int_p)lda_,
	    (float*)x_,
	    (wls_int_p)incx_,
	    (float*)beta_,
	    y_,
	    (wls_int_p)incy_);

    };
  
    template <> inline void blas_t::gemv<double>(const char * trans_,
						 const_int_p m_,
						 const_int_p n_,
						 const double* alpha_,
						 const double*a_,
						 const_int_p lda_,
						 const double*x_,
						 const_int_p incx_,
						 const double* beta_,
						 double*y_,
						 const_int_p incy_)
    {

      dgemv((char*)trans_,
	    (wls_int_p)m_,
	    (wls_int_p)n_,
	    (double*)alpha_,
	    (double*)a_,
	    (wls_int_p)lda_,
	    (double*)x_,
	    (wls_int_p)incx_,
	    (double*)beta_,
	    y_,
	    (wls_int_p)incy_);

    };
  

  };

};
  
