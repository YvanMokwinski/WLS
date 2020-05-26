#pragma once

#include "wls.hpp"
#include "wls-sparse-matrix.hpp"
#include "wls-dense-matrix.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace WLS
{
  void partial_sums(integer_t size_,
		  integer_t * count_);
  
  namespace Input
  {
    //!
    //! @brief Reader for the matrix market format.
    //!
    struct matrix_market_t
    {
    private:      
      FILE * m_in{};
      
    public:
      //!
      //! @brief Get the dimension (number of rows, number of cols.
      //! @param [out] nrows_ 
      //! @param [out] ncols_
      //! @return the status.
      wls_status_t dimension	(wls_int_p 	nrows_,
				 wls_int_p 	ncols_);
      //!
      //! @brief Get the dimension (number of rows, number of cols and number of non-zero elements.
      //! @param [out] nrows_ 
      //! @param [out] ncols_
      //! @param [out] nnz_
      //! @return the status.
      wls_status_t dimension	(wls_int_p 	nrows_,
				 wls_int_p 	ncols_,
				 wls_int_p 	nnz_);
      
      //!
      //! @brief Compute the number of non-zero coefficients following the direction_.
      //! @param storage_ 	
      //! @param nrows_ 	
      //! @param ncols_ 	
      //! @param nnz_ 	
      //! @param count_
      //! @return the status.
      wls_status_t nnz		(mat_storage_t 	storage_,
				 wls_int_t 	nrows_,
				 wls_int_t 	ncols_,
				 wls_int_t 	nnz_,
				 wls_int_p 	count_);

      template<typename T>
      wls_status_t create_dense	(wls_int_t  	nrows_,
				 wls_int_t  	ncols_,
				 T * 		x_,
				 wls_int_t 	ldx_);
      
      template<typename T>
      wls_status_t create_compressed_sparse(mat_storage_t 	storage_,
					    mat_indexing_t	indexing_,
					    wls_int_t  		nrows_,
					    wls_int_t  		ncols_,
					    wls_int_t  		nnz_,
					    wls_int_p 		begin_,
					    wls_int_p 		idx_,
					    T * 		values_);

      ~matrix_market_t();
      static wls_status_t create(matrix_market_t*self_ ,const char * filename_);
      static wls_status_t import(dense::matrix& a_,const char * filename_);
      template<typename T>
      static wls_status_t import(sparse::matrix_t<T>*a_,const char * filename_);
    };
  };
  
};


namespace WLS
{
  namespace Output
  {
    struct matrix_market_t
    {
    protected:
      std::string m_filename;
      std::ofstream m_out;
    public : 
      matrix_market_t(const std::string& filename_);   
      virtual ~matrix_market_t();      
      template <typename _type> matrix_market_t& operator<<(const _type &type_);      
    };
    
    std::ostream& operator << (std::ostream &out_,
			       const WLS::dense::matrix&that_);
  
  };
};






namespace WLS
{
void partial_sums(integer_t size_,
		    integer_t * count_)
  {
    for (integer_t i = 1;i <= size_; ++i)
      {
	count_[i] += count_[i-1];
      }
  };
  namespace Input
  {

    matrix_market_t::~matrix_market_t()
    {
      if (this->m_in)
	{
	  fclose(this->m_in);
	  this->m_in = nullptr;
	}
    };
    
    wls_status_t matrix_market_t::create(matrix_market_t*self_ ,const char * filename_)
    {
      self_->m_in  = fopen(filename_,"r");
      if (!self_->m_in)
	{
	  return status_t::invalid_file;
	}
      return status_t::success;		    
    };

    wls_status_t matrix_market_t::import(dense::matrix& a_,const char * filename_)
      {
	wls_status_t status;
	Input::matrix_market_t matrix_market;
	status = Input::matrix_market_t::create(&matrix_market,filename_);
	if (status_t::success != status)
	  {	    
	    //fprintf(stderr,"status(=" iformat ")", (wls_int_t)status);
	    //	    return status;
	  }
	wls_int_t mat_nrows,mat_ncols;
	status = matrix_market.dimension(&mat_nrows,
					 &mat_ncols);
	if (status_t::success != status)
	  {
	    //      fprintf(stderr,"status(=" iformat ")", (wls_int_t)status);
	    //      return status;
	  }
	double * x = (double*)malloc(sizeof(double)*mat_nrows*mat_ncols);
	dense::matrix::define(a_,mat_nrows,mat_ncols,x,mat_nrows);
	
	status = matrix_market.create_dense(mat_nrows,
					    mat_ncols,
					    x,
					    mat_nrows);
	return status;
      };


    template<typename T>
    wls_status_t matrix_market_t::import(sparse::matrix_t<T>*a_,const char * filename_)
    {
      Input::matrix_market_t matrix_market;
      wls_status_t status;
      
      status = Input::matrix_market_t::create(&matrix_market,filename_);
      if (status_t::success != status)
	{	    
	  //fprintf(stderr,"status(=" iformat ")", (wls_int_t)status);
	  //	    return status;
	}
      wls_int_t mat_nrows,mat_ncols,mat_nnz;
      status = matrix_market.dimension(&mat_nrows,
				       &mat_ncols,
					 &mat_nnz);
	if (status_t::success != status)
	  {
	    //      fprintf(stderr,"status(=" iformat ")", (wls_int_t)status);
	    //      return status;
	  }

	//
	// Allocate.
	//
	wls_int_t*mat_begin 	= (wls_int_t*)calloc((mat_nrows+1),sizeof(wls_int_t));
	if (nullptr == mat_begin)
	  {
	    status = status_t::error_memory;
	    //      goto state_error;
	  }
   
	status = matrix_market.nnz(mat_storage_t::row,
				   mat_nrows,
				   mat_ncols,
				   mat_nnz,
				   mat_begin);
	if (status)
	  {
	    //      fprintf(stderr,"error(=" iformat ")", (wls_int_t)status);
	    //      return status;
	  }
    
	partial_sums(mat_nrows,
			  mat_begin);


	wls_int_t*mat_idx	= (wls_int_t*)malloc(sizeof(wls_int_t)*(mat_nnz));
	if (nullptr == mat_idx)
	  {
	    status = status_t::error_memory;
	    //      goto state_error;
	  }
    
	T*mat_values 	= (T*)malloc(sizeof(T)*(mat_nnz));
	if (nullptr == mat_values)
	  {
	    status = status_t::error_memory;
	    //      goto state_error;
	  }
	
	status = matrix_market.create_compressed_sparse(mat_storage_t::row,
							mat_indexing_t::Fortran,
							mat_nrows,
							mat_ncols,
							mat_nnz,
							mat_begin,
							mat_idx,
							mat_values);
	if (status)
	  {
	    //      std::cerr << "matrix_market_read_cs failed, error code " << status << std::endl;
	    //      return status;
	  }

	a_->init(new sparse::symbolic_t(true,// mat_indexing_t::C==mat_indexing_t::Fortran,
					     mat_nrows,
					     mat_ncols,
					     mat_nnz,
					     mat_begin,
					     mat_idx,
					     true),
		 mat_values);
	return status;
      };
      

      

    wls_status_t matrix_market_t::dimension(wls_int_p 		nrows_,
					wls_int_p 		ncols_,
					wls_int_p 		nnz_)
    {
      if (nullptr == nrows_) 	return status_t::invalid_pointer;
      if (nullptr == ncols_) 	return status_t::invalid_pointer;
      if (nullptr == nnz_) 	return status_t::invalid_pointer;
      
      wls_status_t status = status_t::success;
      
      nrows_[0] = 0;
      ncols_[0] = 0;
      nnz_[0] = 0;
      
      //
      // Read dimensions.
      //
      { char *line = nullptr;
	size_t len;
	for (ssize_t read = getline(&line, &len, m_in);
	     read != -1;
	     read = getline(&line, &len, m_in))
	  {
	    if (line[0] != '%')
	      {
		if (3 != wls_sscanf(line,nrows_,ncols_,nnz_))
		  {
		    status = status_t::invalid_format;
		  }
		break;
	      }
	  } if (line) { free(line); } }
      
      return status;
    };


    
    wls_status_t matrix_market_t::dimension(wls_int_p 		nrows_,
					wls_int_p 		ncols_)
    {
      if (nullptr == nrows_) 	return status_t::invalid_pointer;
      if (nullptr == ncols_) 	return status_t::invalid_pointer;
      
      wls_status_t status = status_t::success;
      
      nrows_[0] = 0;
      ncols_[0] = 0;
      
      //
      // Read dimensions.
      //
      { char *line = nullptr;
	size_t len;
	for (ssize_t read = getline(&line, &len, m_in);
	     read != -1;
	     read = getline(&line, &len, m_in))
	  {
	    if (line[0] != '%')
	      {
		if (2 != wls_sscanf(line,nrows_,ncols_))
		  {
		    status = status_t::invalid_format;
		  }
		break;
	      }
	  } if (line) { free(line); } }
      
      return status;
    };

  //!
  //! @brief Compute the number of non-zero coefficients following the direction_.
  //! @param in_ The file descriptor.
  //! @param direction_ 	The direction.
  //! @param nrows_ 		The direction.
  //!
  wls_status_t matrix_market_t::nnz	(mat_storage_t 	storage_,
				 wls_int_t 	nrows_,
				 wls_int_t 	ncols_,
				 wls_int_t 	nnz_,
				 wls_int_p 	count_)
  {
    fpos_t pos;
    fgetpos(m_in,&pos);

    if (nullptr == count_) 	return status_t::invalid_pointer;
    if (nrows_ < 0) 		return status_t::invalid_size;
    if (ncols_ < 0) 		return status_t::invalid_size;
    if (nnz_ < 0)   		return status_t::invalid_size;
    wls_status_t 	status 		= status_t::success;
    wls_int_t 	nnz 		= 0;  
    { size_t len;
      char *line = nullptr;
      wls_int_t ij[2];
      for (ssize_t read = getline(&line, &len, m_in);
	   read != -1;
	   read = getline(&line, &len, m_in))
	{
	  if (line[0] != '%')
	    {
	      if (2 != wls_sscanf(line, &ij[0], &ij[1]))
		{
		  status = status_t::invalid_format;
		  break;
		}
	      else if (ij[0] <= 0 || ij[0] > nrows_ || ij[1] <= 0 || ij[1] > ncols_)
		{
#ifndef NDEBUG
		  fprintf(stderr,"// diagonostic, error from invalid index line %d of source file\n",__LINE__);
#endif
		  status = status_t::invalid_index;
		  break;
		}
	      else
		{
		  count_[ij[storage_]-1+1] += 1;
		  ++nnz;
		}
	    }
	} if (line) { free(line); } }
    
    if (status_t::success != status)
      {
	return status;
      }

    if (nnz != nnz_)
      {
	return status_t::invalid_size;
      }

    fsetpos(m_in,&pos);
    return status;  
  };

  
    template<typename T>
    wls_status_t matrix_market_t::create_dense(wls_int_t  		nrows_,
					       wls_int_t  		ncols_,
					       T * 			x_,
					       wls_int_t 		ldx_)
    {

      if (nullptr == x_) 	return status_t::invalid_pointer;
      if (nrows_ < 0) 		return status_t::invalid_size;
      if (ncols_ < 0) 		return status_t::invalid_size;
      if (ldx_ < nrows_) return status_t::invalid_size;
      { size_t len;
	char * line = nullptr;
	for (wls_int_t irow=0;irow<nrows_;++irow)
	  {
	    ssize_t read = getline(&line, &len, this->m_in);
	    if (read==-1)
	      {
		return status_t::invalid_format;
	      }
	    for (wls_int_t icol=0;icol<ncols_;++icol)
	      {
		if (1 != wls_sscanf(line,&x_[ldx_ * icol + irow]))
		  {
		    return status_t::invalid_format;
		  }
	      }  
	  } if (line) free(line); }
      return status_t::success;
    };

    template<typename T>
    wls_status_t matrix_market_t::create_compressed_sparse(mat_storage_t 	storage_,
							   mat_indexing_t	indexing_,
							   wls_int_t  		nrows_,
							   wls_int_t  		ncols_,
							   wls_int_t  		nnz_,
							   wls_int_p 		begin_,
							   wls_int_p 		idx_,
							   T * 			values_)
    {

      if (nullptr == begin_) 	return status_t::invalid_pointer;
      if (nullptr == idx_) 	return status_t::invalid_pointer;
      if (nullptr == values_) 	return status_t::invalid_pointer;
      if (nrows_ < 0) 		return status_t::invalid_size;
      if (ncols_ < 0) 		return status_t::invalid_size;
      if (nnz_ < 0)   		return status_t::invalid_size;
      if (mat_storage_t::is_invalid(storage_)) return status_t::invalid_enum;
    
      wls_int_t idx_base = indexing_.base();
      wls_status_t  status   = status_t::success;
      wls_int_t nrowcols = (mat_storage_t::row == storage_) ? nrows_ : ncols_;
    
      { size_t len;
	char * line = nullptr;
	wls_int_t indices[2];
	T value;
	for (ssize_t read = getline(&line, &len, this->m_in);
	     read != -1;
	     read = getline(&line, &len, this->m_in))
	  {
	    if (line[0] != '%')
	      {	  
		if (3 != wls_sscanf(line,&indices[0],&indices[1], &value))
		  {
		    status = status_t::invalid_format;
		    break;
		  }
		else if (indices[0] <= 0 || indices[0] > nrows_ || indices[1] <= 0 || indices[1] > ncols_)
		  {
		    fprintf(stderr,"" iformat " " iformat " " iformat " " iformat "\n",indices[0],nrows_,indices[1],ncols_);
		    status = status_t::invalid_index;
		    break;
		  }
		else
		  {
		    wls_int_t at    	= begin_[indices[storage_] - 1]++;		
		    idx_[at] 		= indices[ (storage_+1)%2 ] + (idx_base - 1);
		    values_[at] 		= value;
		  }
	      }
	  } if (line) { free(line); } }
  
      if (status_t::success!=status)
	{
	  return status;
	}
  
      for (wls_int_t i = nrowcols-1;i>0;--i)
	{
	  begin_[i] = begin_[i-1];
	}
      begin_[0] = 0;
  
      {
    
	//
	// Now we sort the coefficients.
	//   
	wls_int_t mx = begin_[1]-begin_[0];
	for (wls_int_t i = 1;i<nrowcols;++i)
	  {
	    wls_int_t n = begin_[i+1]-begin_[i];
	    mx = (mx < n) ? n : mx;
	  }
    
	void * mem = malloc(2*sizeof(double)*mx);
#ifndef NDEBUG
	std::cout << "mx " << mx << " / " << nrowcols << std::endl;
#endif
	for (wls_int_t i = 0;i<nrowcols;++i)
	  {
	    wls_int_t bound = begin_[i+1];
	    for (wls_int_t j = begin_[i];j<bound;++j)
	      {
		T * m  	= (T*)mem;
		T* x 	= (T*)(m + 2* (j-begin_[i]) + 1);
		wls_int_t* k = (wls_int_t*)(x - 1);
		*k = idx_[j];
		*x = values_[j];
	      }
	
	    qsort(mem,begin_[i+1]-begin_[i],2*sizeof(double),
		  [](const void * a,const void * b)
		  {
		    wls_int_t*a_ = (wls_int_t*)a;
		    wls_int_t*b_ = (wls_int_t*)b;
		    if (a_[0] > b_[0])
		      {
			return 1;
		      }
		    else if (a_[0] < b_[0])
		      {
			return -1;
		      }
		    else
		      {
			return 0;
		      }
		  });
	
	    for (wls_int_t j = begin_[i];j<begin_[i+1] ;++j)
	      {
		T * m  = (T*)mem;
		T* x = (T*)(m + 2* (j-begin_[i]) + 1);
		wls_int_t* k = (wls_int_t*)(x - 1);
		idx_[j] = *k;
		values_[j] = *x;
	      }
	
	  }
	free(mem);
      }
  
      //
      // begin is in 'C'-indexed.
      //
      if (mat_indexing_t::Fortran == indexing_)
	{
	  for (wls_int_t i = 0;i<=nrowcols;++i)
	    {
	      begin_[i]+=1;
	    }
	}

#if 0      
      for (wls_int_t i = 0;i< nrowcols;++i)
	{
	  for (wls_int_t j = begin_[i];j< begin_[i+1];++j)
	    {
	      printf("%Ld %Ld %Ld\n",i+1,idx_[j-1],j);
	    }
	  printf("\n");
	}
#endif
      //
      // Normal exit.
      //
      return status;

    };

  };

};




namespace WLS
{
  namespace Output
  {

    matrix_market_t::matrix_market_t(const std::string& filename_)
      : m_filename(filename_),
	m_out(m_filename.c_str()) 
    {
      this->m_out.precision(15);
      this->m_out.setf(std::ios::scientific);
    };
    
    matrix_market_t::~matrix_market_t()
    {
      this->m_out.close();
    };
    
    template <typename _type> matrix_market_t& matrix_market_t::operator<<(const _type &type_)
    {
      m_out << type_;
      return *this;
    };
        
    std::ostream& operator << (std::ostream &out_,
			       const dense::matrix&that_)    
    {
      out_ << that_.m << " " << that_.n << std::endl;
      for (wls_int_t i=0;i < that_.m;++i)
	{
	  for (wls_int_t j=0;j<that_.n;++j)
	    {
	      out_ << " " << that_.x[that_.ld*j+i];
	    }
	  out_ << std::endl;
	}    
      return out_;
    };
  
  };
};
