#pragma once

#include "WLS/WLSConfig.hpp"
#include <string> 
namespace WLS
{

  //! @brief Abstract base class for all solvers.
  //! An InverseOperator computes the solution of \f$ A(x)=b\f$ where                                    
  //! \f$ A : X \to Y \f$ is an operator.
  //! Note that the solver "knows" which operator
  //! to invert and which preconditioner to apply (if any). The
  //! user is only interested in inverting the operator.
  //! InverseOperator might be a Newton scheme, a Krylov subspace method, or a direct solver or just anything.
  class IInverseOperator
  {
  public:
  
    //! 
    //! @brief The error message if any operation failed.
    //! 
    virtual std::string GetErrorMessage() const = 0;
  
    //! 
    //! @brief The size of the required temporary array of double to apply the linear operator.
    //! 
    //! @return The size of the required temporary vector to apply the linear operator.</returns>
    virtual WLS::integer_t GetSizeOfTemporaryVector() const = 0;
  
    //! 
    //! @brief Compute the inverse operator.
    //! 
    //! @param outHasFailed Indicates if computing the inverse operator has failed. 
    //! True if it has failed, false otherwise.
    virtual void Compute(bool* outHasFailed) = 0;

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
    virtual void Apply(const char*transpose,
		       double*__restrict__ y,
		       const double * __restrict__ x,
		       const WLS::integer_t  tmpSize,
		       double*__restrict__  tmp,
		       bool* outHasFailed) = 0;
    
  };

};
