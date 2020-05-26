#pragma once
#include "wls-sparse-blas.hpp"
#include "wls-mkl.hpp"

namespace WLS
{
  namespace sparse
  {
    
    template <> void blas_t::csrgemv1based<double>(const char* transa_,
						   const blas_t::integer_t m_,
						   const double * a_,
						   const blas_t::integer_pt ia_,
						   const blas_t::integer_pt ja_,
						   const double*x_,
						   double*y_)
    {

      mkl_dcsrgemv((char*)transa_,
		   (blas_t::integer_pt)&m_,
		   (double*)a_,
		   (blas_t::integer_pt)ia_,
		   (blas_t::integer_pt)ja_,
		   (double*)x_,
		   (double*)y_);

    };

    template <> void blas_t::csrgemv1based<float>(const char* transa_,
						  const blas_t::integer_t m_,
						  const float * a_,
						  const blas_t::integer_pt ia_,
						  const blas_t::integer_pt ja_,
						  const float*x_,
						  float*y_)
    {

      mkl_scsrgemv((char*)transa_,
		   (blas_t::integer_pt)&m_,
		   (float*)a_,
		   (blas_t::integer_pt)ia_,
		   (blas_t::integer_pt)ja_,
		   (float*)x_,
		   (float*)y_);

    };



    template <> void blas_t::csrtrsv0based<double>(const char* uplo_,
						   const char* transa_,
						   const char* diag_,
						   const blas_t::integer_t m_,
						   const double* a_,
						   const blas_t::integer_pt ia_,
						   const blas_t::integer_pt ja_,
						   const double* x_,
						   double* y_)
    {

      mkl_cspblas_dcsrtrsv((char*)uplo_,
			   (char*)transa_,
			   (char*)diag_,
			   (blas_t::integer_pt)&m_,
			   (double*)a_,
			   (blas_t::integer_pt) ia_,
			   (blas_t::integer_pt) ja_,
			   (double*) x_,
			   (double*) y_);

    };
  
    template <> void blas_t::csrtrsv0based<float>(const char* uplo_,
						  const char* transa_,
						  const char* diag_,
						  const blas_t::integer_t m_,
						  const float* a_,
						  const blas_t::integer_pt ia_,
						  const blas_t::integer_pt ja_,
						  const float* x_,
						  float* y_)
    {

      mkl_cspblas_scsrtrsv((char*)uplo_,
			   (char*)transa_,
			   (char*)diag_,
			   (blas_t::integer_pt)&m_,
			   (float*)a_,
			   (blas_t::integer_pt) ia_,
			   (blas_t::integer_pt) ja_,
			   (float*) x_,
			   (float*) y_);

    };


    template <> void blas_t::csrtrsv1based<double>(const char* uplo_,
						   const char* transa_,
						   const char* diag_,
						   const blas_t::integer_t m_,
						   const double* a_,
						   const blas_t::integer_pt ia_,
						   const blas_t::integer_pt ja_,
						   const double* x_,
						   double* y_)
    {

      mkl_dcsrtrsv((char*)uplo_,
		   (char*)transa_,
		   (char*)diag_,
		   (blas_t::integer_pt)&m_,
		   (double*)a_,
		   (blas_t::integer_pt) ia_,
		   (blas_t::integer_pt) ja_,
		   (double*) x_,
		   (double*) y_);

    };

    template <> void blas_t::csrtrsv1based<float>(const char* uplo_,
						  const char* transa_,
						  const char* diag_,
						  const blas_t::integer_t m_,
						  const float* a_,
						  const blas_t::integer_pt ia_,
						  const blas_t::integer_pt ja_,
						  const float* x_,
						  float* y_)
    {

      mkl_scsrtrsv((char*)uplo_,
		   (char*)transa_,
		   (char*)diag_,
		   (blas_t::integer_pt)&m_,
		   (float*)a_,
		   (blas_t::integer_pt) ia_,
		   (blas_t::integer_pt) ja_,
		   (float*) x_,
		   (float*) y_);

    };

  };

};
