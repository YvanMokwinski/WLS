#pragma once

#include <array>
#include <iostream>

#include "wls.h"

#ifdef WLS_WITH_MKL

#include "wls-dense-blas-mkl.hpp"
#include "wls-dense-lapack-mkl.hpp"

#else

#endif

namespace WLS
{
  
#ifdef WLS_ILP64
  using integer_t  = long long int;
#else
  using integer_t  = int;
#endif
  using integer_pt = integer_t*__restrict__;

  struct mat_indexing_t
  {
  public: typedef enum ekind : integer_t
    { C = 0, Fortran } kind;
    
  private: kind m_value;
    
  public: inline  mat_indexing_t(integer_t k)
    : m_value((kind)k)
    {
    };
    
  public: inline  mat_indexing_t()
    : m_value(C)
    {
    };
    
  public: inline operator integer_t() const
    {
      return this->m_value;
    };
    
  public: inline integer_t base() const
    {
      return (mat_indexing_t::Fortran == m_value) ? 1 : 0;
    } 
    
  };

  
  struct mat_storage_t
  {
  public: typedef enum ekind : integer_t
    { row = 0, col } kind;
    
  private: kind m_value;
    
  public: inline  mat_storage_t(integer_t k)
    : m_value((kind)k)
    {
    };
    
  public: inline  mat_storage_t()
    : m_value(row)
    {
    };
    
  public: inline operator integer_t() const
    {
      return this->m_value;
    };
  public: inline static bool is_invalid(mat_storage_t value_) 
    {
      return value_ != row && value_ != col;
    };
  };

  struct status_t
  {
  public: typedef enum ekind : integer_t
    { 	success = 0,
	error_memory,
	invalid_index,
	invalid_file,
	invalid_format,
	invalid_pointer,
	invalid_size,
	invalid_enum } kind;
    
  private: kind m_value;
    
  public: inline  status_t(integer_t k)
    : m_value((kind)k)
    {
    };
    
  public: inline  status_t()
    : m_value(success)
    {
    };
    
  public: inline operator integer_t() const
    {
      return this->m_value;
    };
    
  };


  //!
  //! @brief Interface for a linear operator.
  //!
  class linear_operator
  {
  public:
    //!    
    //! @brief Apply the linear operator.
    //! @param transpose_ Transpose.
    //! @param y_ The output vector.
    //! @param x_ The input vector.
    //!
    virtual void Apply(const char * 			transpose_,
		       double*__restrict__ 		y_,
		       const double*__restrict__ 	x_) = 0;
  };


  //! @brief Abstract base class for all solvers.
  //! An InverseOperator computes the solution of \f$ x=A^{-1}b\f$ by solving the linear system \f$ A x =b\f$ where 
  //! \f$ A : X \to Y \f$ is an operator.
  //! Note that the solver "knows" which operator
  //! to invert and which preconditioner to apply (if any). 
  //! The user is only interested in inverting the operator.
  //! InverseOperator might be a Newton scheme, a Krylov subspace method, or a direct solver or just anything.
  class inverse_operator
  {
  public:
  
    //! 
    //! @brief The error message if any operation failed.
    //! 
    virtual std::string get_error_message() const = 0;
  
    //! 
    //! @brief The size of the required temporary array of double to apply the linear operator.
    //! 
    //! @return The size of the required temporary vector to apply the linear operator.</returns>
    virtual WLS::integer_t get_buffer_size() const = 0;
  
    //! 
    //! @brief Compute the inverse operator.
    //! 
    //! @param outHasFailed Indicates if computing the inverse operator has failed. 
    //! True if it has failed, false otherwise.
    virtual void compute(bool* outHasFailed) = 0;

    //! 
    //! Apply the inverse operator: y := Op(x).
    //! 
    //! @param transpose No transpose ('N' or 'n'), Transpose ('T' or 't').
    //! @param y The output vector.
    //! @param x The input vector.
    //! @param tmpSize The size of the temporary vector.
    //! @param tmp The temporary array of double.
    //! @param outHasFailed Indicates if applying the operator has failed. 
    //! True if it has failed, false otherwise.
    virtual void apply(const char*transpose,
		       double*__restrict__ y,
		       const double * __restrict__ x,
		       const WLS::integer_t  tmpSize,
		       double*__restrict__  tmp,
		       bool* outHasFailed) = 0;
    
  };

  
};

  inline int wls_sscanf(const char * str_,WLS::integer_t * i_,WLS::integer_t * j_)
  {
    return sscanf(str_,iformat " " iformat, i_, j_);
  };

  inline int wls_sscanf(const char * str_,WLS::integer_t * m_,WLS::integer_t * n_,WLS::integer_t * nnz_)
  {
    return sscanf(str_,iformat " " iformat " " iformat, m_, n_, nnz_);
  };

  template<typename real_t>
  inline int wls_sscanf(const char * str_,WLS::integer_t * i_,WLS::integer_t * j_,real_t*x_);

  template<typename real_t>
  inline int wls_sscanf(const char * str_,real_t*x_);

  template<>
  inline int wls_sscanf<double>(const char * str_,double*x_)
  {
    return sscanf(str_, "%le" , x_);
  }

  template<>
  inline int wls_sscanf<float>(const char * str_,float*x_)
  {
    return sscanf(str_, "%e" , x_);
  }

  template<>
  inline int wls_sscanf<double>(const char * str_,WLS::integer_t * i_,WLS::integer_t * j_,double*x_)
  {
    return sscanf(str_,iformat " " iformat " %le" , i_, j_, x_);
  }

  template<>
  inline int wls_sscanf<float>(const char * str_,WLS::integer_t * i_,WLS::integer_t * j_,float*x_)
  {
    return sscanf(str_,iformat " " iformat " %e" , i_, j_, x_);
  }

