#pragma once

#include "WLSConfig.hpp"
#include "WLA/include/matrix.hpp"

namespace WLS
{
  
    
  //! 
  //! @brief Implementation of a diagonal block matrix with homogeneous blocks.
  //! 
  class DiagonalBlockMatrix
  {
    
    //! 
    //! @brief The number of blocks of the diagonal matrix.
    //! 
  public: WLS::integer_t m_nblocks;
    
    //! 
    //! @brief The size of a block diagonal matrix.
    //! 
  public: WLS::integer_t m_size_block;
    
    //! 
    //! @brief The size of the diagonal matrix.
    //! 
  protected: WLS::integer_t m_size;
    
    //!
    //! @brief The block matrices.
    //!
  public: WLA::matrix m_matrix;
    
    //!
    //! @brief Constructor.
    //! @param nblocks_ The number of diagonal blocks.
    //! @param size_block_ The size of a diagonal block.
    //!
  public: DiagonalBlockMatrix(const WLS::integer_t nblocks_,
			      const WLS::integer_t size_block_)
    : m_nblocks(nblocks_),
      m_size_block(size_block_),
      m_size(nblocks_ * size_block_),
      m_matrix(size_block_ * size_block_,
	       nblocks_)
    {
    };
    
    //!
    //! Destructor.
    //! 
  public: virtual ~DiagonalBlockMatrix()
    {
    };
    
    //!
    //! @brief Get/Set the value of the diagonal coefficient.
    //!
    //! @param index The index of the diagonal coefficient.
    //! @return The value of the diagonal coefficient.</returns>
    //!
  private: const double* __restrict__ operator[](WLS::integer_t block_index_) const 
    {
      return this->m_matrix(0,block_index_);
    };
    
  private: double* __restrict__ operator[](WLS::integer_t block_index_) 
    {
      return this->m_matrix(0,block_index_);
    };
    
  private: static void block_amux(const char* 		transpose_,
				  WLA::matrix_h& 	block_,
				  double 		rx_,
				  WLA::matrix_h& 	x_,
				  double 		ry_,
				  WLA::matrix_h& 	y_)
    {
      if (transpose_[0] == 't' || transpose_[0] == 'T')
	{
	  if (1.0 == ry_)
	    {
	      y_ += block_.transpose() * x_;
	    }
	  else if (0.0 == ry_)
	    {
	      y_ = block_.transpose() * x_;
	    }
	  else
	    {
	      std::cerr<< "not implemented" << std::endl;
	      exit(1);
	    }
	}
      else
	{
	  if (1.0 == ry_)
	    {
	      y_ += rx_ * block_ * x_;
	    }	  
	  else if (0.0 == ry_)
	    {
	      y_ = rx_ * block_ * x_;
	    }
	  else
	    {
	      std::cerr<< "not implemented" << std::endl;
	      exit(1);
	    }
	}
    };

  public: void operator() (WLS::integer_t i,WLA::matrix_h& h)
    {
      h.m = this->m_size_block;
      h.n = h.m;
      h.ld = h.m;
      h.x = static_cast<WLA::matrix_h&>(this->m_matrix)(0,i);
    };
    
    //!
    //! Perform the matrix-vector product.
    //!
    //! @param transpose_ No transpose ('N' or 'n'), Transpose ('T' or 't').
    //! @param  y_ The output matrix.
    //! @param  x_ The input matrix.
    //!
  public: void amux(const char* transpose_,
		    WLA::matrix_h& x_,
		    WLA::matrix_h& y_)
    {
      WLA::matrix_h block;      
      WLA::matrix_h x = x_.subcols(0,1);      
      WLA::matrix_h y = y_.subcols(0,1);      
      for (WLS::integer_t block_index = 0; block_index < this->m_nblocks; ++block_index)
	{
	  WLA::matrix_h::define(block,m_size_block,m_size_block,(*this)[block_index],m_size_block);
	  x.x = x_(0, block_index);
	  y.x = y_(0, block_index);
	  block_amux(transpose_,
		     block,
		     1.0,
		     x,
		     0.0,
		     y);
	}
    };

    //!
    //! Perform the matrix-vector product.
    //!
    //! @param transpose_ No transpose ('N' or 'n'), Transpose ('T' or 't').
    //! @param  y_ The output vector.
    //! @param  x_ The input vector.
    void Amux(const char* transpose_,
	      double*__restrict__ 	y_,
	      const double*__restrict__ x_)
    {
      
      WLA::matrix_h y(this->m_size_block,this->m_nblocks,y_,this->m_size_block);
      WLA::matrix_h x(this->m_size_block,this->m_nblocks,(double*)x_,this->m_size_block);
      this->amux(transpose_,x,y);
    };
    
    //!
    //! @brief Compute the inverse.
    //!
    //! @param minAbsoluteValue Threshold to apply on the absolute value of the diagonal coefficient,
    //! the default value is MathAndStats.MachinePrecisionDouble.
    //!
    void block_inverse(double * __restrict__ workspace_,
		       const double minAbsoluteValue = ((double)1.1e-16))
    {
      const WLS::integer_t n 	= this->m_size_block;

      WLA::matrix_h W(n,n,workspace_,n);
      WLA::matrix_h D(n,n,(*this)[0],n);
      
      WLA::matrix Id = WLA::matrix::identity(n);

      WLS::integer_t * lcperm = new WLS::integer_t[n];
      for (WLS::integer_t block_index=0;block_index<this->m_nblocks;++block_index)
	{
	  D.x 	= (*this)[block_index];
	  W 	= D;
	  D     = Id;
	  
	  WLS::integer_t info_lapack;      
	  dgesv(&n,
		&n,
		W.x,
		&W.ld,
		lcperm,
		D.x,
		&D.ld,
		&info_lapack);
	  if (info_lapack != 0)
	    {
	      //
	      // Factorization failed, let's replace it.
	      //
	      std::cerr << "WLS::DiagonalBlockMatrix::Inverse failed, block " << block_index << " " << std::endl;
	      D = Id;
	    }
	}      
      delete[]lcperm;
    };

#if 0     


#endif    

  };

  std::ostream& operator << (std::ostream&out_,const DiagonalBlockMatrix&that_)
  {
    out_ << "Diagonal matrix " << std::endl;
    out_ << that_.m_matrix << std::endl;
    return out_;
  };
  
  //! 
  //! @brief Implementation of a diagonal block matrix with homogeneous blocks.
  //! 
  class DiagonalBlockMatrixFactory
  {
    
  private: DiagonalBlockMatrix * m_diagonal_matrix;    
  public:  DiagonalBlockMatrix * create_inverse()
    {
      //      std::cout << *m_diagonal_matrix << std::endl;
      WLS::integer_t n = m_diagonal_matrix->m_size_block;
      WLA::matrix_h D;
      
      WLA::matrix WW 		= WLA::matrix(n,n);
      WLA::matrix_h& W 		= static_cast<WLA::matrix_h&>(WW);
      
      WLS::integer_t info_lapack;
      WLS::integer_t * lcperm 	= new WLS::integer_t[n];
      WLA::matrix id 		= WLA::matrix::identity(n);
      for (WLS::integer_t block_index=0;block_index< m_diagonal_matrix->m_nblocks;++block_index)
	{	  
	  (*m_diagonal_matrix)(block_index, D);
	  //	  static_cast<WLA::matrix_h&>(W) = D;
	  W = D;
	  D = id;	  
	  dgesv(&n,
		&n,
		W.x,
		&W.ld,
		lcperm,
		D.x,
		&D.ld,
		&info_lapack);
	  
	  if (info_lapack != 0)
	    {
	      //
	      // Factorization failed, let's replace it.
	      //
	      std::cerr << "WLS::DiagonalBlockMatrix::Inverse failed, block " << block_index << " " << std::endl;
	      exit(1);
	      //	      D = Id;
	    }
	}
      delete [] lcperm;
      auto p = this->m_diagonal_matrix;
      this->m_diagonal_matrix = nullptr;
      return p;
    };
    
  public: DiagonalBlockMatrixFactory() = delete;
  public: DiagonalBlockMatrixFactory(const WLS::integer_t nblocks_,
				     const WLS::integer_t size_block_)
    {
      this->m_diagonal_matrix = new DiagonalBlockMatrix(nblocks_,
							size_block_); 
    };

  public: void addelm(WLS::integer_t 	blocki_,
		      double          	s_,
		      WLA::matrix_h&  	block_matrix_h_)
    {
      
      WLA::matrix_h D;

      (*m_diagonal_matrix)(blocki_,D);
      D = s_ * block_matrix_h_;
    };
    
  public: void addelm(WLS::integer_t 	blocki_,
		      double          	s_,
		      double *  	blockx_,
		      WLS::integer_t 	blockld_)
    {
      WLS::integer_t n = m_diagonal_matrix->m_size_block;
      WLA::matrix_h block(n,n,blockx_,blockld_);
      addelm(blocki_,s_,block);      
    };
  };
  
};
