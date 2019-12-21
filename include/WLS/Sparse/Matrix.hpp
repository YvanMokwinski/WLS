#pragma once

#include "WLS/WLSConfig.hpp"
#include "WLS/ILinearOperator.hpp"
#include "WLS/Sparse/SparsityPattern.hpp"
#include "WLS/Sparse/OpMatrixAssembly.hpp"
#include "WLS/Sparse/SpBlas.hpp"

namespace WLS
{
  namespace Sparse
  {
    template <typename _real_t> class Matrix
    {

    private: using real_t = _real_t;
    private: using this_t = Matrix<real_t>;
    
    public: Matrix(const SparsityPattern * sparsityPattern_);

    public: inline bool 			HasFortranIndexing() const noexcept;
    public: inline void 			clear() noexcept;    
    public: inline WLS::integer_t 		GetN() 	const noexcept;
    public: inline WLS::integer_t 		GetNC() const noexcept;
    public: inline WLS::integer_pt 		GetB() 	const noexcept;  
    public: inline WLS::integer_pt  		GetI() 	const noexcept;

    public: inline real_t*__restrict__ 		GetX() noexcept;
    public: inline const real_t*__restrict__ 	GetX() const noexcept;
    public: inline const SparsityPattern * GetSparsityPattern() const noexcept
      {
	return this->m_sparsityPattern;
      };
      
    public: inline Matrix& operator = (const real_t value_) noexcept;
      
    public: inline bool Ass(const WLS::integer_t	i_,
			    const WLS::integer_t	j_,
			    const real_t	x_) noexcept;
  
    public: template <typename _indirect_array_t>
    inline OpMatrixAssembly<this_t,_indirect_array_t> operator[] (const _indirect_array_t& indirection_) noexcept;

    private: const SparsityPattern * m_sparsityPattern;
      
    private: real_t * __restrict__ m_x; // WLS::Vector<real_t> m_x;

    public:
      class LinearOperator : public WLS::ILinearOperator
      {
      private:
	const this_t& m_matrix;
      public:
	LinearOperator(const this_t&matrix_)
	  : m_matrix(matrix_)
	{	
	};
      
      public:
	virtual void Apply(const char * transpose_,
			   double*__restrict__ y_,
			   const double*__restrict__ x_)
	{
	  WLS::Sparse::SpBlas::csrgemv1based(transpose_,
					     this->m_matrix.GetN(),
					     this->m_matrix.GetX(),
					     this->m_matrix.GetB(),
					     this->m_matrix.GetI(),
					     x_,
					     y_);	
	};
      };
    
    public:
      WLS::ILinearOperator* GetLinearOperator() const
      {
	return new LinearOperator(*this);
      };
    
    
    };

    template <typename real_t>
    inline bool Matrix<real_t>::HasFortranIndexing() const noexcept
      {
	return this->m_sparsityPattern->HasFortranIndexing();
      };

    
    template <typename real_t> Matrix<real_t>::Matrix(const SparsityPattern * sparsityPattern_)
      : m_sparsityPattern(sparsityPattern_) // ,
	//	m_x(m_sparsityPattern->GetNC())    
    {
      m_x = (real_t*)malloc(sizeof(real_t) * m_sparsityPattern->GetNC());
    };
  
    template <typename real_t>
    template<typename _indirect_array_t>
    inline OpMatrixAssembly<Matrix<real_t>,_indirect_array_t> Matrix<real_t>::operator[] (const _indirect_array_t& indirection_) noexcept
    {
      return {*this,indirection_};
    };
  
    template <typename real_t> inline Matrix<real_t>& Matrix<real_t>::operator = (const real_t value_) noexcept
    {
      auto nc = m_sparsityPattern->GetNC();
      for (integer_t i = 0;i<nc;++i)
	{
	  this->m_x[i] = value_;
	}
      
      //      m_x = value_;
    };
  
    template <typename real_t> inline bool Matrix<real_t>::Ass(const WLS::integer_t	i_,
							       const WLS::integer_t	j_,
							       const real_t		x_) noexcept
    {
      const WLS::integer_t where = this->m_sparsityPattern->Find(i_,j_);
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
  
    template <typename real_t> inline void Matrix<real_t>::clear() noexcept
    {
      auto nc = m_sparsityPattern->GetNC();
      for (integer_t i = 0;i<nc;++i)
	{
	  this->m_x[i] = 0.0;
	}
    };
  
    template <typename real_t> inline WLS::integer_t Matrix<real_t>::GetN() const noexcept
    {
      return this->m_sparsityPattern->GetN();
    };
  
    template <typename real_t> inline WLS::integer_t Matrix<real_t>::GetNC() const noexcept
    {
      return this->m_sparsityPattern->GetNC();
    };
  
    template <typename real_t> inline WLS::integer_pt Matrix<real_t>::GetB() const noexcept
    {
      return this->m_sparsityPattern->GetB();
    };
  
    template <typename real_t> inline WLS::integer_pt  Matrix<real_t>::GetI() const noexcept
    {
      return this->m_sparsityPattern->GetI();
    };
  
    template <typename real_t> inline real_t*__restrict__ Matrix<real_t>::GetX() noexcept
    {
      return this->m_x;
    };

    template <typename real_t> inline const real_t*__restrict__ Matrix<real_t>::GetX() const noexcept
    {
      return this->m_x;
    };

  };
};
