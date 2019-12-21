#pragma once

#include "WLSConfig.hpp"

namespace WLS
{
  //! 
  //! @brief Implementation of a diagonal matrix.
  //! 
  class DiagonalMatrix
  {
  protected :    
    //! 
    //! @brief The size of the diagonal matrix.
    //! 
    WLS::integer_t m_size;
    
    //!
    //! @brief The values of the diagonal matrix.
    //!
    double*__restrict__ m_values;    
  public:
    
    //!
    //! @brief Constructor.
    //! @param size The size of the diagonal matrix.
    //!
    DiagonalMatrix(const WLS::integer_t size_)
    {
      this->m_size = size_;
      this->m_values = new double[size_];
    };
    
    //!
    //! Destructor.
    //! 
    virtual ~DiagonalMatrix()
    {
      if (nullptr != this->m_values)
	{
	  delete [] this->m_values;
	  this->m_values = nullptr;
	}
      this->m_size = 0;
     };
    
    //!
    //! @brief Get/Set the value of the diagonal coefficient.
    //!
    //! @param index The index of the diagonal coefficient.
    //! @return The value of the diagonal coefficient.</returns>
    double operator[](const WLS::integer_t index) const 
    {
      return m_values[index];
    };
    
    //!
    //! @brief Compute the inverse.
    //!
    //! @param minAbsoluteValue Threshold to apply on the absolute value of the diagonal coefficient,
    //! the default value is MathAndStats.MachinePrecisionDouble.
    //!
    void Inverse(const double minAbsoluteValue = ((double)1.1e-16))
    {
      static constexpr const double s_zero = double(0.0);
      static constexpr const double s_one = double(1.0);
      for (WLS::integer_t index = 0; index < this->m_size; ++index)
	{
	  double v = this->m_values[index];	   
	  if (((v < s_zero ? -v : v) < minAbsoluteValue))
	    {
	      v = (v < s_zero)
		? -minAbsoluteValue
		: minAbsoluteValue;
	    }	  
	  this->m_values[index] = s_one / v;
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
      if ( (transpose_[0] == 'N') ||
	   (transpose_[0] == 'n') ||
	   (transpose_[0] == 'T') ||
	   (transpose_[0] == 't'))
	{
	  for (WLS::integer_t index = 0; index < this->m_size; ++index)
	    {
	      y_[index] = x_[index] * this->m_values[index];
	    }
	}   
    };
    
  };
  
};
