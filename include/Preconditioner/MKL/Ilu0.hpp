#pragma once

#include "WLS/ILinearOperator.hpp"
#include "WLS/Iterative/MKL/Parameters.hpp"

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_rci.h"

namespace WLS
{
  namespace Preconditioner
  {
    namespace MKL
    {
  
      /// <summary>
      /// Implementation of the Incomplete LU Factorization with level 0 of fill-in.
      /// </summary>
      /// <remarks>
      /// Avoid using this preconditioner with the Conjugate Gradient sparse iterative solver because in general, 
      /// it produces a non-symmetric resulting matrix even if the original matrix is symmetric.
      /// </remarks>
      class Ilu0 : public IInverseOperator
      {
      private:

	struct Error
	{
	public:
	  /// <summary>
	  /// Enumerate error values.
	  /// </summary>
	  typedef enum _Value
	    {
	      /// <summary>
	      /// Indicates that the task completed normally.
	      /// </summary>
	      Completed = 0,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted and that error occurred: 
	      /// at least one diagonal element is omitted from the matrix in CSR3 format 
	      /// </summary>
	      MissingDiagonalElement = -101,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted because the matrix contains a diagonal element with the value of zero.
	      /// </summary>
	      ZeroDiagonalElement = -102,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted because the matrix contains a diagonal element which is so small that it could cause an overflow, 
	      /// or that it would cause a bad approximation to ILU0.
	      /// </summary>
	      TooSmallDiagonalElement = -103,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted because the memory is insufficient for the internal work array.
	      /// </summary>
	      NotEnoughInternalMemory = -104,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted because the size n of the input matrix is less than 0.
	      /// </summary>
	      InvalidMatrixDimension = -105,
		
	      /// <summary>
	      /// Indicates that the routine was interrupted because the column indices ja are not in the ascending order.
	      /// </summary>
	      InvalidMatrixIndexing = -106,
	    } Value;

	public:
	  /// <summary>
	  /// Get the error message from the routine <see cref="ILUT"/>.
	  /// </summary>
	  /// <param name="error">The parameter 'ierr' from the routine.</param>
	  /// <returns>The error message related to the error.</returns>
	  static std::string GetMessage(const Error::Value error)
	  {
	    switch (error)
	      {
	      case Error::Completed:
		{
		  return
		    "Indicates that the task completed normally.";
		}

	      case Error::MissingDiagonalElement:
		{
		  return
		    "Indicates that the routine was interrupted and that error occurred: at least one diagonal element is omitted from the matrix in CSR3 format.";
		}
		 
	      case Error::ZeroDiagonalElement:
		{
		  return "Indicates that the routine was interrupted because the matrix contains a diagonal element with the value of zero.";
		}
		 
	      case Error::TooSmallDiagonalElement:
		{
		  return
		    "Indicates that the routine was interrupted because the matrix contains a diagonal element which is so small that it could cause an overflow, or that it would cause a bad approximation to ILU0.";
		}
		 
	      case Error::NotEnoughInternalMemory:
		{
		  return
		    "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
		}

	      case Error::InvalidMatrixDimension:
		{
		  return
		    "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
		}

	      case Error::InvalidMatrixIndexing:
		{
		  return
		    "Indicates that the routine was interrupted because the column indices ja are not in the ascending order.";
		}

	      default:
		{
		  return "Unknown error value";
		}
	      }
	  };
	};
	
      private:
	/// <summary>
	/// The matrix A.
	/// </summary>
	Sparse::Matrix<double> * m_A;
	  
	/// <summary>
	/// The preconditioner.
	/// </summary>
	double * __restrict__ m_ilu0Values;
	  
	/// <summary>
	/// Parameter from RCI FGMRES, first initialized by the FGMRES method.
	/// </summary>
	WLS::Iterative::MKL::Parameters * m_parameters;
	  
	/// <summary>
	/// Status.
	/// </summary>
	Error::Value m_error;
	  
      public:    

	virtual ~Ilu0()
	{
	  this->m_A = NULL;
	  if (NULL != this->m_ilu0Values)
	    {
	      free(this->m_ilu0Values);
	      this->m_ilu0Values = NULL;
	    }
	  this->m_parameters = NULL;
	  this->m_error = Error::Completed;	  
	};

	/// <summary>
	/// Create an Incomplete LU Factorization preconditioner with a level-0 of fill-in.
	/// </summary>
	/// <param name="A">The matrix A from which we need to compute a preconditioner.</param>
	/// <param name="parameters">Parameters of the RCI FGMRES computations, first initialized by the MKL FGMRES.</param>
	Ilu0(Sparse::Matrix<double>* A,
	     WLS::Iterative::MKL::Parameters* parameters)
	{
	    
	  //
	  // Assign data members
	  //
	  this->m_parameters 	= parameters;
	  this->m_A 		= A;
	  this->m_ilu0Values 	= (double * __restrict__)malloc(sizeof(double)*A->GetNC());

	  WLS::integer_t* paramIntegers 	= parameters->GetParamIntegers();
	  double * __restrict__ paramReals 	= parameters->GetParamReals();
	    
	  //
	  // ipar[30]
	  // specifies how the routine operates when a zero diagonal element occurs during calculation. 
	  // If this parameter is set to 0 (the default value set by the routine dfgmres_init), then the calculations are stopped 
	  // and the routine returns a non-zero error value. Otherwise, the diagonal element is set to the value of dpar[31] and the calculations continue.
	  // 
	  paramIntegers[30] = 1;
	    
	  //
	  // specifies a small value, which is compared with the computed diagonal elements. 
	  // When ipar[30] is not 0, then diagonal elements less than dpar[30] are set to dpar[31]. 
	  // The default value is 1.0e-16.
	  //   
	  // Note
	  // This parameter can be set to the negative value, because the calculation uses its absolute value.
	  // If this parameter is set to 0, the comparison with the diagonal element is not performed.
	  //
	  paramReals[30] = 1.0e-16;
	    
	  //
	  // specifies the value that is assigned to the diagonal element if its value is less 
	  // than dpar[30] (see above). The default value is 1.0e-10.
	  //
	  paramReals[31] = 1.0e-10;
	    
	};
    
      public:
	  
	/// <see cref="IInverseOperator.SizeOfTemporaryVector"/>>
	WLS::integer_t GetSizeOfTemporaryVector()const
	{
	  return this->m_A->GetN();
	};

	/// <see cref="IInverseOperator.ErrorMessage"/>
	std::string GetErrorMessage()const
	{
	  return Error::GetMessage(m_error); 
	};
	  
	/// <see cref="IInverseOperator.Compute"/>
	void Compute(bool* outHasFailed)
	{
	  *outHasFailed = false;
	  this->m_error = Error::Completed;
	    
	  //
	  // Set the matrix to 1-based indexing.
	  //
	  const bool hasFortranIndexing = m_A->HasFortranIndexing();
	  if (false == hasFortranIndexing)
	    {
	      std::cerr << "need to come here" << std::endl;
	      exit(1);
	      //	      Sparse_fortran_indexation(this->m_A);
	    }
	    
	  //
	  // Compute the incomplete factorization
	  //
	    
	  {
	    WLS::integer_t n = this->m_A->GetN();
	    WLS::integer_t lierr;
	    dcsrilu0(&n,
		     this->m_A->GetX(),
		     this->m_A->GetB(),
		     this->m_A->GetI(),
		     this->m_ilu0Values,
		     m_parameters->GetParamIntegers(),
		     m_parameters->GetParamReals(),
		     &lierr);
	    this->m_error = (Error::Value)lierr;
	  }
	    
	  //
	  // Re-set the matrix to 0-based indexing.
	  //
	  if (false == hasFortranIndexing)
	    {	      
	      // Sparse_c_indexation(this->m_A);
	    }
	    
	  switch (this->m_error)
	    {
	    case Error::Completed:
	      {
		break;
	      }
	    default:
	      {
		*outHasFailed = true;
		break;
	      }
	    }
	    
	};
	  
	  
	/// <see cref="IInverseOperator.Apply"/>
	void Apply(const char*	transpose,
		   double * __restrict__ 		y,
		   const double * __restrict__		x,
		   const WLS::integer_t  		tmpSize,
		   double * __restrict__  		tmp,
		   bool* 		outHasFailed)
	{
      
	  *outHasFailed = false;
	  WLS::integer_t dimension = this->m_A->GetN();

	  const char * sL = "L";
	  const char * sU = "U";
	  const char * sN = "N";
	  const bool hasFortranIndexing = m_A->HasFortranIndexing();
	  if (hasFortranIndexing)
	    {
	      //
	      // Forward substitution
	      //
	      mkl_dcsrtrsv(sL,
			   sN,
			   sU,
			   &dimension,
			   this->m_ilu0Values,
			   this->m_A->GetB(),
			   this->m_A->GetI(),
			   x,
			   tmp);
		
		
	      //
	      // Backward substitution
	      //
	      mkl_dcsrtrsv(sU,
			   sN,
			   sN,
			   &dimension,
			   this->m_ilu0Values,
			   m_A->GetB(),
			   m_A->GetI(),			   
			   tmp,
			   y);
	    }
	  else
	    {
		
	      //
	      // Forward substitution
	      //
	      mkl_cspblas_dcsrtrsv(sL,
				   sN,
				   sU,
				   &dimension,
				   this->m_ilu0Values,
				   m_A->GetB(),
				   m_A->GetI(),			   
				   x,
				   tmp);
		
		
	      //
	      // Backward substitution
	      //
	      mkl_cspblas_dcsrtrsv(sU,
				   sN,
				   sN,
				   &dimension,
				   this->m_ilu0Values,
				   m_A->GetB(),
				   m_A->GetI(),			   
				   tmp,
				   y);
	    }
	    
	};
	  
	  
      };

    }; // namespace MKL

  }; // namespace Preconditioner
  
}; // namespace WLS
