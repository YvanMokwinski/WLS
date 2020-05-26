#pragma once

#include "wla_traits.hpp"
#include "wla_cache_op.hpp"

namespace WLS
{
  namespace dense
  {
  //! 
  //! @brief Struct to define a matrix handle with a column-oriented storage.
  //! @remark
  //! - A matrix handle is meta-data information on how accessing
  //! the data of a matrix given its pointer to the first element, the number of rows,
  //! the number of columns and the leading dimension.
  //! - A matrix handle is not responsible of the memory.
  //!
  struct matrix_h
  {
    double*__restrict__  x{};
    wla_int_t   m{};
    wla_int_t   n{};
    wla_int_t   ld{};
    
    
    //!
    //! @brief Empty constructor.
    //!
    matrix_h();

    //!
    //! @brief Destructor.
    //!
    ~matrix_h();

    
    //!
    //! @brief Constructor.
    //! @param m_ The number of rows.
    //! @param n_ The number of columns.
    //! @param x_ The pointer to the first element.
    //! @param ld_ The leaging dimension.
    //!
    matrix_h(wla_int_t m_, wla_int_t n_,double*__restrict__  x_,wla_int_t ld_);

    static void define(matrix_h&that_,wla_int_t m_, wla_int_t n_,double*__restrict__  x_,wla_int_t ld_) 
    {
      assert(m_ > 0);
      assert(n_ > 0);
      assert(nullptr != x_);
      assert(ld_ >= m_);
      that_.x=x_;
      that_.m=m_;
      that_.n=n_;
      that_.ld = ld_;
    };
    
    
    //!
    //! @brief Get the pointer to the element (i_,j_).
    //! @param i_ The row index.
    //! @param j_ The column index.
    //! @return The pointer to the specified element.
    //!
    inline double * __restrict__ operator ()(wla_int_t i_,wla_int_t j_);
    
    //!
    //! @brief Get the constant pointer to the element (i_,j_).
    //! @param i_ The row index.
    //! @param j_ The column index.
    //! @return The constant pointer to the specified element.
    //!
    inline const double * __restrict__ operator ()(wla_int_t i_,wla_int_t j_) const;
    
    inline matrix_h subh(wla_int_t srow_,
			 wla_int_t erow_,
			 wla_int_t scol_,
			 wla_int_t ecol_);
    
    inline matrix_h subrows(wla_int_t srow_,
			    wla_int_t erow_);
    
    inline matrix_h subcols(wla_int_t scol_,
			    wla_int_t ecol_);

    
    inline void lus(matrix_h&rhs_,wla_int_t*__restrict__ ipiv)
    {
      wla_int_t info_lapack;
      lapack_t::getrs("N",
			&this->m,
			&rhs_.n,
			this->x,
			&this->ld,
			ipiv,
			rhs_(0,0),
			&rhs_.ld,
			&info_lapack);
    };
    
    inline void lu(wla_int_t*__restrict__ ipiv)
    {
      wla_int_t info_lapack;
      lapack_t::getrf(&this->m,
		      &this->m,
		      this->x,
		      &this->ld,
		      ipiv,
		      &info_lapack);      
    };

    //  rhs->x,
    //  &h->ld,
    
    //!
    //! @brief Assigned addition operator.
    //! @param that_ That matrix handle.
    //! @return This.
    //!
    inline matrix_h& operator = (const matrix_h& that_);

    //!
    //! @brief Assigned addition operator.
    //! @param that_ That matrix handle.
    //! @return This.
    //!
    inline matrix_h& operator += (const matrix_h& that_);
    //!
    //! @brief Assigned substract operator.
    //! @param that_ That matrix handle.
    //! @return This.
    //!
    inline matrix_h& operator -= (const matrix_h& that_);


    inline matrix_h& operator =  (double alpha_);
    inline matrix_h& operator *= (double alpha_);
    
    inline void clear();
    
    void gesv(struct vector_h& rhs,wla_int_t * __restrict__ lcperm);
    
    inline temp_transpose transpose() const;    
    inline matrix_h& operator+= (const temp_gemm&temp);
    inline matrix_h& operator-= (const temp_gemm&temp);
    inline matrix_h& operator = (const temp_gemm&temp);
    inline matrix_h& operator = (const temp_scal&temp);
    inline temp_gemv operator * (const vector_h&v_) const;  
    inline temp_gemm operator * (const matrix_h&v_) const;  
    inline temp_scal operator * (const double&v_) const;
  };
  
  temp_scal operator * (const double&v_,const matrix_h &m) 
  {
    return {v_,'N',m};
  };

  inline matrix_h matrix_h::subh(wla_int_t srow_,
				 wla_int_t erow_,
				 wla_int_t scol_,
				 wla_int_t ecol_) 
  {
    
    assert(srow_ + 1 > 0 && srow_ <= this->m);
    assert(erow_ > srow_ && erow_ <= this->m);
    
    assert(scol_ + 1 > 0 && scol_ <= this->n);
    assert(ecol_ > scol_ && ecol_ <= this->n);
    
    return matrix_h(erow_-srow_,
		    ecol_-scol_,
		    (*this)(srow_,scol_),
		    this->ld);
  };
  
  inline matrix_h matrix_h::subrows(wla_int_t srow_,
				    wla_int_t erow_) 
  {
    assert(srow_ + 1 > 0 && srow_ <= this->m);
    assert(erow_ > srow_ && erow_ <= this->m);
    return matrix_h(erow_-srow_,
		    this->n,
		    (*this)(srow_,0),
		    this->ld);
  };
  
  inline matrix_h matrix_h::subcols(wla_int_t scol_,
				    wla_int_t ecol_) 
  {
    assert(scol_ + 1 > 0 && scol_ <= this->n);
    assert(ecol_ > scol_ && ecol_ <= this->n);
    return matrix_h(this->m,
		    ecol_-scol_,
		    (*this)(0,scol_),
		    this->ld);
  };
  
  //!
  //! @brief Copy operator.
  //! @param that_ That matrix handle.
  //! @return This.
  //!
  inline matrix_h& matrix_h::operator = (const matrix_h& that_)
  {
    static constexpr wla_int_t n1 = 1;    
    assert(that_.n == this->n);
    assert(that_.m == this->m);      
    if ( (this->ld == this->m) && (this->ld == that_.ld) )
      {
	wla_int_t nn = that_.n*that_.ld;
	 blas_t::copy(&nn, that_(0,0), &n1, (*this)(0,0),&n1);
      }
    else
      {
	wla_int_t nn = (that_.m < this->m) ? that_.m : this->m;
	for (wla_int_t j=0;j<this->n;++j)
	  {
	    blas_t::copy(&nn, that_(0,j), &n1, (*this)(0,j),&n1);
	  }	    
      }
    return *this;
  };

  //!
  //! @brief Assigned addition operator.
  //! @param that_ That matrix handle.
  //! @return This.
    //!
    inline matrix_h& matrix_h::operator += (const matrix_h& that_)
    {
      static constexpr double r1    = 1.0;
      static constexpr wla_int_t n1 = 1;    
      assert(that_.n == this->n);
      assert(that_.m == this->m);      
      if ( (this->ld == this->m) && (this->ld == that_.ld) )
	{
	  wla_int_t nn = that_.n*that_.ld;
	  blas_t::axpy(&nn,&r1, that_(0,0), &n1, (*this)(0,0),&n1);
	}
      else
	{
	  wla_int_t nn = (that_.m < this->m) ? that_.m : this->m;
	  for (wla_int_t j=0;j<this->n;++j)
	    {
	      blas_t::axpy(&nn,&r1, that_(0,j), &n1, (*this)(0,j),&n1);
	    }	    
	}
      return *this;
    };

    //!
    //! @brief Assigned substract operator.
    //! @param that_ That matrix handle.
    //! @return This.
    //!
    inline matrix_h& matrix_h::operator -= (const matrix_h& that_)
    {
      static constexpr double mr1    = -1.0;
      static constexpr wla_int_t n1 = 1;    
      assert(that_.n == this->n);
      assert(that_.m == this->m);      
      if ( (this->ld == this->m) && (this->ld == that_.ld) )
	{
	  wla_int_t nn = that_.n*that_.ld;
	  blas_t::axpy(&nn,&mr1, that_(0,0), &n1, (*this)(0,0),&n1);
	}
      else
	{
	  wla_int_t nn = (that_.m < this->m) ? that_.m : this->m;
	  for (wla_int_t j=0;j<this->n;++j)
	    {
	      blas_t::axpy(&nn,&mr1, that_(0,j), &n1, (*this)(0,j),&n1);
	    }	    
	}
      return *this;
    };

    inline matrix_h& matrix_h::operator = (double alpha_)
    {
      for (wla_int_t j=0;j<this->n;++j)
	{
	  for (wla_int_t i=0;i<this->m;++i)
	    {    
	      (*this)(i,j) [0] = alpha_;
	    }
	}
      return *this;
    };
    
    
  inline matrix_h& matrix_h::operator *= (double alpha_)
    {
      wla_int_t n1 = 1;
      if (this->m == ld)
	{
	  wla_int_t N = this->m*this->n;
	  blas_t::scal(&N,&alpha_, x, &n1);
	}
      else
	{
	  for (wla_int_t j=0;j<this->n;++j)
	    {
	      blas_t::scal(&this->m,&alpha_, (*this)(0, j), &n1);
	    }
	}
      return *this;
    };

  inline temp_transpose matrix_h::transpose() const
    {
      return {'T',*this};
    };
    
  inline void matrix_h::clear()
    {
      if (ld == this->m)
	{
	  const wla_int_t N = this->n * this->m;
	  for (wla_int_t i=0;i<N;++i)
	    {
	      x[i] = 0.0;
	    }
	}
      else
	{
	  for (wla_int_t j=0;j<this->n;++j)
	    {
	      for (wla_int_t i=0;i<this->m;++i)
		{    
		  x[j*ld+i] = 0.0;
		}
	    }
	}
    }

  
  inline matrix_h& matrix_h::operator=(const temp_scal&temp)
  {
    *this = temp.A;
    if (temp.trans == 'T')
      {
	/* transpose first */
	if (this->n == this->m)
	  {
	    wla_int_t N = this->n;
	    for (wla_int_t j=0;j<N;++j)
	      {
		for (wla_int_t i=j+1;i<N;++i)
		  {
		    auto tmp = x[j*ld+i];
		    x[j*ld+i] = x[i*ld+j];
		    x[i*ld+j] = tmp;
		  }
	      }
	  }
	else
	  {
	    // 0 1 2  3 4 5  6 7 8  9 10 11
	    // 0 3 6 9   1 4 7 10   2 5 8 11
	    // x0 x1 x2 x3;
	    // x1;
	    // x2;
	    std::cerr << "not implemented" << std::endl;
	    exit(1);
	  }	
      }
    *this *= temp.a;
    return *this;
  };

  inline matrix_h& matrix_h::operator+=(const temp_gemm&temp)
  {
#ifndef NDEBUG
      wla_int_t a_n = (temp.transA == 'N') ? temp.A.n : temp.A.m;
      wla_int_t a_m = (temp.transA == 'N') ? temp.A.m : temp.A.n;
      wla_int_t b_n = (temp.transB == 'N') ? temp.B.n : temp.B.m;
      wla_int_t b_m = (temp.transB == 'N') ? temp.B.m : temp.B.n;
      assert(a_n==b_m);
      assert(a_m==this->m);
      assert(b_n==this->n);
#endif
      static constexpr double r1 = 1.0;
      wla_int_t k = (temp.transA == 'N') ? temp.A.n : temp.A.m;
      blas_t::gemm(&temp.transA,
		    &temp.transB,
		    &this->m,
		    &this->n,
		    &k,
		    &temp.a,
		    temp.A.x,
		    &temp.A.ld,
		    temp.B.x,
		    &temp.B.ld,
		    &r1,
		    x,
		    &ld);
      
      return *this;
    };
    
  inline matrix_h& matrix_h::operator-=(const temp_gemm&temp)
    {
#ifndef NDEBUG
      wla_int_t a_n = (temp.transA == 'N') ? temp.A.n : temp.A.m;
      wla_int_t a_m = (temp.transA == 'N') ? temp.A.m : temp.A.n;
      wla_int_t b_n = (temp.transB == 'N') ? temp.B.n : temp.B.m;
      wla_int_t b_m = (temp.transB == 'N') ? temp.B.m : temp.B.n;
      assert(a_n==b_m);
      assert(a_m==this->m);
      assert(b_n==this->n);
#endif
      static constexpr double r1 = 1.0;
      double a = -temp.a;     
      wla_int_t k = (temp.transA == 'N') ? temp.A.n : temp.A.m;
      blas_t::gemm(&temp.transA,
		    &temp.transB,
		    &this->m,
		    &this->n,
		    &k,
		    &a,
		    temp.A.x,
		    &temp.A.ld,
		    temp.B.x,
		    &temp.B.ld,
		    &r1,
		    x,
		    &ld);
      return *this;
    };
    


  
  inline matrix_h& matrix_h::operator=(const temp_gemm&temp)
  {
#ifndef NDEBUG
    wla_int_t a_n = (temp.transA == 'N') ? temp.A.n : temp.A.m;
    wla_int_t a_m = (temp.transA == 'N') ? temp.A.m : temp.A.n;
    wla_int_t b_n = (temp.transB == 'N') ? temp.B.n : temp.B.m;
    wla_int_t b_m = (temp.transB == 'N') ? temp.B.m : temp.B.n;
    assert(a_n==b_m);
    assert(a_m==this->m);
    assert(b_n==this->n);
#endif
    static constexpr double r0 = 0.0;
    wla_int_t k = (temp.transA == 'N') ? temp.A.n : temp.A.m;
    blas_t::gemm(&temp.transA,
		  &temp.transB,
		  &this->m,
		  &this->n,
		  &k,
		  &temp.a,
		  temp.A.x,
		  &temp.A.ld,
		  temp.B.x,
		  &temp.B.ld,
		  &r0,
		  x,
		  &ld);
    return *this;
  };
  
  matrix_h::matrix_h(){};
  matrix_h::~matrix_h()
  {
  };

  
  double * __restrict__ matrix_h::operator ()(wla_int_t i_,wla_int_t j_)
  {
    return this->x + (j_*this->ld + i_);
  };
  
  const double * __restrict__ matrix_h::operator ()(wla_int_t i_,wla_int_t j_) const
  {
    return this->x + (j_*this->ld + i_);
  };
  
  matrix_h::matrix_h(wla_int_t m_, wla_int_t n_,double*__restrict__  x_,wla_int_t ld_) : x(x_), m(m_), n(n_), ld(ld_)
  {
    assert(m_ > 0);
    assert(n_ > 0);
    assert(nullptr != x_);
    assert(ld_ >= m_);
  };
  
  inline temp_gemv matrix_h::operator * (const vector_h&v_) const
  {
    return {1.0,'N',*this, v_};
  };
  
  inline temp_gemm matrix_h::operator * (const matrix_h&v_) const
  {
    return {1.0,'N','N',*this, v_};
  };

  inline temp_scal matrix_h::operator * (const double&v_) const
  {
    return {v_,'N',*this};
  };
  
  std::ostream& operator << (std::ostream&out_,const matrix_h&that_)
  {
    out_ << "[ ";
    wla_int_t mm1 = that_.m-1;
    for (wla_int_t i=0;i < mm1;++i)
      {
	if (i > 0)
	  {
	    out_ << "  ";
	  }
	out_ << that_.x[that_.ld*0+i];
	for (wla_int_t j=1;j<that_.n;++j)
	  {
	    out_ << ", " << that_.x[that_.ld*j+i];
	  }
	out_ << ";" << std::endl;
      }    
    out_ << "  " << that_.x[that_.ld*0+mm1];
    for (wla_int_t j=1;j<that_.n;++j)
      {
	out_ << ", " << that_.x[that_.ld*j+mm1];
      }
    out_ << " ];";
    return out_;
  };



  struct vector_h : public matrix_h
  {
    vector_h() : matrix_h()
    {

    };
    vector_h(wla_int_t n_, double*__restrict__  x_,wla_int_t incx_)
      : matrix_h(1,n_,x_,incx_)
    {};

    
    static void define(vector_h&that_,wla_int_t n_,double*__restrict__  x_,wla_int_t incx_) 
    {
      assert(n_ > 0);
      assert(nullptr != x_);
      assert(incx_ >= 1);
      that_.x=x_;
      that_.m=1;
      that_.n=n_;
      that_.ld = incx_;
    };

    template <typename F>
    inline void apply(F f)
    {
      const wla_int_t n = this->n;
      for (wla_int_t j=0;j<n;++j)
	{
	  double*__restrict__ e = x + j*ld;
	  *e = f(*e);
	}    
    };

    inline vector_h& operator=(const temp_gemv&temp)
    {
      static constexpr double r0 = 0.0;
#ifndef NDEBUG
      wla_int_t a_n = (temp.trans == 'N') ? temp.A.n : temp.A.m;
      wla_int_t a_m = (temp.trans == 'N') ? temp.A.m : temp.A.n;
      wla_int_t x_n     = temp.b.n;
      assert(a_n==x_n);
      assert(a_m==this->n);
#endif
      
      blas_t::gemv(&temp.trans,
		    &temp.A.m,
		    &temp.A.n,
		    &temp.a,
		    temp.A.x,
		    &temp.A.ld,
		    temp.b.x,
		    &temp.b.ld,
		    &r0,
		    this->x,
		    &this->ld);
      return *this;
    };
  
    inline vector_h& operator+=(const temp_gemv&temp)
    {
#ifndef NDEBUG
      wla_int_t a_n = (temp.trans == 'N') ? temp.A.n : temp.A.m;
      wla_int_t a_m = (temp.trans == 'N') ? temp.A.m : temp.A.n;
      wla_int_t x_n     = temp.b.n;
      assert(a_n==x_n);
      assert(a_m==this->n);
#endif
      static constexpr double r1 = 1.0;
      blas_t::gemv(&temp.trans,
		    &temp.A.m,
		    &temp.A.n,
		    &temp.a,
		    temp.A.x,
		    &temp.A.ld,
		    temp.b.x,
		    &temp.b.ld,
		    &r1,
		    this->x,
		    &this->ld);
      return *this;

    };

    
    inline vector_h& operator-=(const temp_gemv&temp)
    {
#ifndef NDEBUG
      wla_int_t a_n = (temp.trans == 'N') ? temp.A.n : temp.A.m;
      wla_int_t a_m = (temp.trans == 'N') ? temp.A.m : temp.A.n;
      wla_int_t x_n = temp.b.n;
      assert(a_n==x_n);
      assert(a_m==this->n);
#endif
      static constexpr double r1 = 1.0;
      double a = -temp.a;

      blas_t::gemv(&temp.trans,
		   &temp.A.m,
		   &temp.A.n,
		   &a,
		   temp.A.x,
		   &temp.A.ld,
		   temp.b.x,
		   &temp.b.ld,
		   &r1,
		   this->x,
		   &this->ld);
      
      return *this;
    };

  };

  void matrix_h::gesv(struct vector_h& rhs,wla_int_t * __restrict__ lcperm)
  {
    wla_int_t n1=1;
    wla_int_t info_lapack;
    lapack_t::gesv(&this->n,
		    &n1,
		    this->x,
		    &this->ld,
		    lcperm,
		    rhs.x,
		    &this->n,
		    &info_lapack);
  };
  
  template <wla_int_t _m,wla_int_t _n>  struct smatrix : public matrix_h
  {
  private:double s_mem[_m*_n];
  public: smatrix() : matrix_h(_m,_n,s_mem,_m) {};
  };
  
  struct matrix : public matrix_h
  {
  public:static matrix identity(wla_int_t m_)
    {
      matrix id(m_,m_);
      static_cast<matrix_h&>(id) = 0.0;
      for (wla_int_t i=0;i<m_;++i)
	{
	  id.x[i*id.ld+i]=1.0;
	}
      return id;
    };
    
  private: double * __restrict__ m_mem{};
  public: operator matrix_h&() {return *this;};
  public: operator const matrix_h&() const {return *this;};
  private: static double*alloc(size_t size_)
    {
      return (double*__restrict__)malloc(sizeof(double)*size_);
    };
  public: matrix(){};
  public: matrix(wla_int_t m_,
		 wla_int_t n_) : matrix_h(m_,n_,alloc(m_*n_),m_)
    {
      this->m_mem = this->x;
      //  m_mem = (double*__restrict__)malloc(sizeof(double)*m_*n_);
      //  this->x = m_mem;
    };

    static void define(matrix&that_,wla_int_t m_, wla_int_t n_,double*__restrict__  x_,wla_int_t ld_) 
    {
      matrix_h::define(that_,m_,n_,x_,ld_);
      that_.m_mem = x_;
    };
    
    ~matrix()
    {
      if (this->m_mem)
	{
	  free(this->m_mem);
	  this->m_mem = nullptr;
	  this->x = nullptr;
	}
    };
    
  };

};

};
