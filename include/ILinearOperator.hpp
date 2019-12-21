#pragma once

namespace WLS
{
  //!
  //! @brief Interface for a linear operator.
  //!
  class ILinearOperator
  {
  public:
    //!    
    //! @brief Apply the linear operator.
    //! @param transpose_ Transpose.
    //! @param y_ The input vector.
    //! @param x_ The output vector.
    //!
    virtual void Apply(const char * transpose_,
		       double*__restrict__ y_,
		       const double*__restrict__ x_) = 0;
  };
  
};
