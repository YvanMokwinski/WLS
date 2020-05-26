#pragma once

#include "wls.hpp"
#include "wls-io.hpp"
#include "wls-sparse-symbolic.hpp"
#include "wls-sparse-assembly.hpp"
#include "wls-sparse-blas-mkl.hpp"

namespace WLS
{
  namespace sparse
  {
    template <typename T> class matrix_t
    {      
    private: using this_t = matrix_t<T>;
      
    public: matrix_t(){};
    public: matrix_t(const symbolic_t * sparsityPattern_);
    public: matrix_t(const symbolic_t * sparsityPattern_,T*x_);
    public: inline bool 			HasFortranIndexing() const noexcept;
    public: inline void 			clear() noexcept;    
    public: inline wls_int_t 		GetN() 	const noexcept;
    public: inline wls_int_t 		GetNC() const noexcept;
    public: inline wls_int_p 		GetB() 	const noexcept;  
    public: inline wls_int_p  		GetI() 	const noexcept;
      
    public: inline T*__restrict__ 		GetX() noexcept;
    public: inline const T*__restrict__ 	GetX() const noexcept;
      
      //!
      //! @brief Get the (read-only) sparsity pattern.
      //!
    public: inline const symbolic_t * GetSparsityPattern() const noexcept
      {
	return this->m_sparsityPattern;
      };
      
      //!
      //! @brief Assign a scalar.
      //!
    public: inline matrix_t& operator = (T value_) noexcept;
    public: void init(const symbolic_t * sparsityPattern_,T*x_)
      {
	m_sparsityPattern=sparsityPattern_;
	m_x=x_;
      };
      
    public: inline bool Ass(const wls_int_t	i_,
			    const wls_int_t	j_,
			    T	x_) noexcept;
  
    public: template <typename _indirect_array_t>
    inline assembly_t<this_t,_indirect_array_t> operator[] (const _indirect_array_t& indirection_) noexcept;

    private: const symbolic_t * m_sparsityPattern;
      
    private: T * __restrict__ m_x; // WLS::Vector<T> m_x;

    public:
      class operator_mvp : public WLS::linear_operator
      {
      private:
	const this_t& m_matrix;
      public:
	operator_mvp(const this_t&matrix_)
	  : m_matrix(matrix_)
	{	
	};
      
      public:
	virtual void Apply(const char * transpose_,
			   double*__restrict__ y_,
			   const double*__restrict__ x_)
	{
	  WLS::sparse::blas_t::csrgemv1based(transpose_,
					     this->m_matrix.GetN(),
					     this->m_matrix.GetX(),
					     this->m_matrix.GetB(),
					     this->m_matrix.GetI(),
					     x_,
					     y_);	
	};
      };
    
    public:
      WLS::linear_operator * GetLinearOperator() const
      {
	return new operator_mvp(*this);
      };
    
    
    };

    template <typename T>
    inline bool matrix_t<T>::HasFortranIndexing() const noexcept
      {
	return this->m_sparsityPattern->HasFortranIndexing();
      };

    
    template <typename T> matrix_t<T>::matrix_t(const symbolic_t * symbolic_,T*x_)
      : m_sparsityPattern(symbolic_),m_x(x_)
    {
    };

    template <typename T> matrix_t<T>::matrix_t(const symbolic_t * symbolic_)
      : m_sparsityPattern(symbolic_)
    {
      m_x = (T*)malloc(sizeof(T) * m_sparsityPattern->GetNC());
    };
  
    template <typename T>
    template<typename _indirect_array_t>
    inline assembly_t<matrix_t<T>,_indirect_array_t> matrix_t<T>::operator[] (const _indirect_array_t& indirection_) noexcept
    {
      return {*this,indirection_};
    };
    
    template <typename T> inline matrix_t<T>& matrix_t<T>::operator = (T value_) noexcept
    {
      auto nc = m_sparsityPattern->GetNC();
      for (integer_t i = 0;i<nc;++i)
	{
	  this->m_x[i] = value_;
	}
    };
  
    template <typename T> inline bool matrix_t<T>::Ass(const wls_int_t	i_,
							       const wls_int_t	j_,
							       T			x_) noexcept
    {
      const wls_int_t where = this->m_sparsityPattern->Find(i_,j_);
      if (where >= 0)
	{
	  this->m_x[where] += x_;
	  return true;
	}
      else
	{
	  return false;
	}
    };
  
    template <typename T> inline void matrix_t<T>::clear() noexcept
    {
      auto nc = m_sparsityPattern->GetNC();
      for (integer_t i = 0;i<nc;++i)
	{
	  this->m_x[i] = 0.0;
	}
    };
  
    template <typename T> inline wls_int_t matrix_t<T>::GetN() const noexcept
    {
      return this->m_sparsityPattern->GetN();
    };
  
    template <typename T> inline wls_int_t matrix_t<T>::GetNC() const noexcept
    {
      return this->m_sparsityPattern->GetNC();
    };
  
    template <typename T> inline wls_int_p matrix_t<T>::GetB() const noexcept
    {
      return this->m_sparsityPattern->GetB();
    };
  
    template <typename T> inline wls_int_p  matrix_t<T>::GetI() const noexcept
    {
      return this->m_sparsityPattern->GetI();
    };
  
    template <typename T> inline T*__restrict__ matrix_t<T>::GetX() noexcept
    {
      return this->m_x;
    };

    template <typename T> inline const T*__restrict__ matrix_t<T>::GetX() const noexcept
    {
      return this->m_x;
    };

  };
};
