#pragma once

#include "wls.hpp"
#include "wls-iterative-mkl-parameters.hpp"

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
      class Ilut : public IInverseOperator
      {

      private:
	struct Error
	{
	  /// <summary>
	  /// Enumerate error values.
	  /// </summary>
	public: typedef enum Value
	  {
	    /// <summary>
	    /// Indicates that the task completed normally.
	    /// </summary>
	    Completed = 0,
		   
	    /// <summary>
	    /// Indicates that the routine was interrupted because of an error: the number of elements 
	    /// in some matrix row specified in the sparse format is equal to or less than 0.
	    /// </summary>
	    EncounteredAnEmptyRow = -101,
		   
	    /// <summary>
	    /// Indicates that the routine was interrupted because the value of the computed diagonal 
	    /// element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar[30]=0.
	    /// </summary>
	    InvalidDiagonalValue = -102,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the element ia[i] is less than or equal to the element ia[i - 1] (see Sparse Matrix Storage Format).
	    /// </summary>
	    EncounteredAnInvalidOffset = -103,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays
	    /// </summary>
	    NotEnoughInternalMemory = -104,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the input value of maxfil is less than 0.
	    /// </summary>
	    InvalidMaxfil = -105,

	    /// <summary>
	    /// Indicates that the routine was interrupted because the size n of the input matrix is less than 0.
	    /// </summary>
	    InvalidMatrixDimension = -106,

	    /// <summary>
	    /// Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n (see Sparse Matrix Storage Format).
	    /// </summary>
	    InvalidMatrixIndexing = -107,

	    /// <summary>
	    /// The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).
	    /// </summary>
	    InvalidMaxfil2 = 101,

	    /// <summary>
	    /// The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)
	    /// </summary>
	    InvalidTolerance = 102,

	    /// <summary>
	    /// The absolute value of tol is greater than value of dpar[30]; it can result in instability of the calculation.
	    /// </summary>
	    IncompatibleToleranceAndDiagonalElementThreshold = 103,

	    /// <summary>
	    /// The value of dpar[30] is equal to 0. It can cause calculations to fail.
	    /// </summary>
	    InvalidDiagonalElementThreshold = 104,
	  } Type;

	  /// <summary>
	  /// Get the error message from the routine <see cref="ILUT"/>.
	  /// </summary>
	  /// <param name="error_">The parameter 'ierr' from the routine.</param>
	  /// <returns>The error message related to the error.</returns>
	public: static std::string GetMessage(const Error::Type error_)
	  {
	    switch (error_)
	      {
	      case Error::Completed:
		{
		  return
		    "Indicates that the task completed normally.";
		}

	      case Error::EncounteredAnEmptyRow:
		{
		  return
		    "Indicates that the routine was interrupted because of an error: the number of elements in some matrix row specified in the sparse format is equal to or less than 0.";
		}

	      case Error::InvalidDiagonalValue:
		{
		  return "Indicates that the routine was interrupted because the value of the computed diagonal element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar[30]=0.";
		}

	      case Error::EncounteredAnInvalidOffset:
		{
		  return
		    "Indicates that the routine was interrupted because the element ia[i] is less than or equal to the element ia[i - 1] (see Sparse Matrix Storage Format)";
		}

	      case Error::NotEnoughInternalMemory:
		{
		  return
		    "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
		}


	      case Error::InvalidMaxfil:
		{
		  return
		    "Indicates that the routine was interrupted because the input value of maxfil is less than 0.";
		}


	      case Error::InvalidMatrixDimension:
		{
		  return
		    "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
		}


	      case Error::InvalidMatrixIndexing:
		{
		  return
		    "Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n (see Sparse Matrix Storage Format).";
		}


	      case Error::InvalidMaxfil2:
		{
		  return
		    "The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).";
		}

	      case Error::InvalidTolerance:
		{
		  return
		    "The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)";
		}

	      case Error::IncompatibleToleranceAndDiagonalElementThreshold:
		{
		  return
		    "The absolute value of tol is greater than value of dpar[30]; it can result in instability of the calculation.";
		}


	      case Error::InvalidDiagonalElementThreshold:
		{
		  return "The value of dpar[30] is equal to 0. It can cause calculations to fail.";
		}


	      default:
		{
		  return "Unknown error value";
		}
	      }
	  };

	};
	
	/// <summary>
	/// The matrix A.
	/// </summary>
      private: Sparse::Matrix<double> *m_A;
	
	/// <summary>
	/// The preconditioner.
	/// </summary>
      private: Sparse::Matrix<double> *m_P;
	
	/// <summary>
	/// Tolerance for threshold criterion for the resulting entries of the preconditioner.
	/// </summary>
      private: double m_tol;
	
	/// <summary>
	/// Maximum fill-in, which is half of the preconditioner bandwidth. 
	/// The number of non-zero elements in the rows of the preconditioner cannot exceed (2*<see cref="m_maxfil"/>+1)
	/// </summary>
      private: wls_int_t m_maxfil;
	
	/// <summary>
	/// Parameter from RCI FGMRES, first initialized by the FGMRES method.
	/// </summary>
      private: WLS::Iterative::MKL::Parameters * m_parameters;
	  
	/// <summary>
	/// Status.
	/// </summary>
      private: Error::Type m_error;
	
	/// <summary>
	/// Get the required of number of coefficients of the incomplete LU factorization.
	/// </summary>
	/// <param name="n_">The dimension.</param>
	/// <param name="maxfil_">The maximum fill.</param>
	/// <returns>The required number of coefficients.</returns>
      private: static wls_int_t GetRequiredNumberOfCoefficients(const wls_int_t n_,
							const wls_int_t maxfil_)
	{
	  return (2 * maxfil_ + 1) * n_ - maxfil_ * (maxfil_ + 1) + 1;
	};

	
	/// <see cref="IInverseOperator.SizeOfTemporaryVector"/>>
      public: wls_int_t GetSizeOfTemporaryVector()const
	{
	  return this->m_A->GetN();
	};

	
	/// <see cref="IInverseOperator.ErrorMessage"/>
      public: std::string GetErrorMessage()const
	{
	  return Error::GetMessage(m_error); 
	};
	
      public: virtual ~Ilut()
	{
	  this->m_A 		= NULL;
	  delete this->m_P;
	  this->m_P = NULL;
	  this->m_parameters 	= NULL;
	  this->m_error 	= Error::Completed;	  
	};

	/// <summary>
	/// Create an Incomplete LU Factorization preconditioner.
	/// </summary>
	/// <param name="A">The matrix A from which we need to compute a preconditioner.</param>
	/// <param name="tol">Tolerance for threshold criterion for the resulting entries of the preconditioner.</param>
	/// <param name="maxfil">Maximum fill-in, which is half of the preconditioner bandwidth. 
	/// The number of non-zero elements in the rows of the preconditioner cannot exceed (2*<paramref name="maxfil"/>+1)</param>
	/// <param name="parameters">Parameters of the RCI FGMRES computations, first initialized by the MKL FGMRES.</param>
      public: Ilut(Sparse::Matrix<double>* A_,
		   const double tol_,
		   const wls_int_t maxfil_,
		   WLS::Iterative::MKL::Parameters * parameters_)
	{
	  
	  //
	  // Assign data members
	  //
	  this->m_tol = tol_;
	  this->m_maxfil = maxfil_;
	  this->m_parameters = parameters_;
	  this->m_A = A_;

	  
	  wls_int_t* paramIntegers 	= this->m_parameters->GetParamIntegers();
	  double * __restrict__ paramReals 	= this->m_parameters->GetParamReals();
	  
	  //
	  // ipar[30]
	  // specifies how the routine operates if the value of the computed diagonal element is less 
	  // than the current matrix row norm multiplied by the value of the parameter tol. 
	  // If ipar[30] = 0, then the calculation is stopped and the routine returns non-zero error value. 
	  // Otherwise, the value of the diagonal element is set to a value determined by dpar[30] 
	  // (see its description below), and the calculations continue.
	  //
	  paramIntegers[30] = 1;
	  
	  //
	  // used to adjust the value of small diagonal elements. 
	  // Diagonal elements with a value less than the current matrix row norm multiplied 
	  // by tol are replaced with the value of dpar[30] multiplied by the matrix row norm.
	  //
	  paramReals[30] = 1.0e-5;

	  //
	  // Get the dimension
	  //
	  const wls_int_t N = this->m_A->GetN();
	  
	  //
	  // Compute the required number of coefficient.
	  //
	  const wls_int_t requiredNumberOfCoefficients = Ilut::GetRequiredNumberOfCoefficients(N,
										       this->m_maxfil);
	  

	  //
	  // The preconditioner P needs to be 1-based indexed.
	  //
	  
	  this->m_P = new Sparse::Matrix<double>(new Sparse::SparsityPattern(true,
								     N,
								     N,
								     requiredNumberOfCoefficients));
#if 0
	  this->m_P = Sparse_new(N,
				 N,
				 requiredNumberOfCoefficients);
	  this->m_P->format = 1;
#endif

	};

	/// <see cref="IInverseOperator.Compute"/>
      public: void Compute(bool* outHasFailed_)
	{
#ifndef WLS_WITH_MKL
	  *outHasFailed_ = true;
#else
	  *outHasFailed_ = false;
	  this->m_error = Error::Completed;
	    
	  //
	  // Set the matrix to 1-based indexing.
	  //
	  const bool hasFortranIndexing = this->m_A->HasFortranIndexing();
	  if (false == hasFortranIndexing)
	    {
	      std::cout << "come here" << std::endl;
	      exit(1);
	      // Sparse_fortran_indexation(this->m_A);
	    }
	    
	  //
	  // Compute the incomplete factorization
	  //
	    
	  {
	    wls_int_t n = this->m_A->GetN();
	    wls_int_t lierr = 0;
#ifdef WLS_WITH_MKL	    
	    dcsrilut(&n,
		     this->m_A->GetX(),
		     this->m_A->GetB(),
		     this->m_A->GetI(),
		     this->m_P->GetX(),
		     this->m_P->GetB(),
		     this->m_P->GetI(),
		     &this->m_tol,
		     &this->m_maxfil,
		     m_parameters->GetParamIntegers(),
		     m_parameters->GetParamReals(),
		     &lierr);
#endif
	    this->m_error = (Error::Type)lierr;
	  }
	    
	  //
	  // Re-set the matrix to 0-based indexing.
	  //
	  if (false == hasFortranIndexing)
	    {
	      std::cout << "come here" << std::endl;
	      exit(1);
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
		*outHasFailed_ = true;
		break;
	      }
	    }
#endif	    
	};

	/// <see cref="IInverseOperator.Apply"/>
      public: void Apply(const char*	transpose_,
			 double * __restrict__		y_,
			 const double * __restrict__	x_,
			 const wls_int_t  	tmpSize_,
			 double * __restrict__ 		tmp_,
			 bool* 		outHasFailed_)
	{
#ifdef WLS_WITH_MKL	          

	  *outHasFailed_ = false;
	  wls_int_t dimension = this->m_A->GetN();

	  static constexpr const char * s_lowerTriangle = "L";
	  static constexpr const char * s_upperTriangle = "U";

	  static constexpr const char * s_notTransposed = "N";
	  
	  static constexpr const char * s_notUnitTriangular = "N";
	  static constexpr const char * s_unitTriangular = "U";


	  //
	  // Forward substitution
	  //
	  mkl_dcsrtrsv(s_lowerTriangle,
		       s_notTransposed,
		       s_unitTriangular,
		       &dimension,
		       this->m_P->GetX(),
		       this->m_P->GetB(),
		       this->m_P->GetI(),
		       x_,
		       tmp_);
	  
	  //
	  // Backward substitution
	  //
	  mkl_dcsrtrsv(s_upperTriangle,
		       s_notTransposed,
		       s_notUnitTriangular,
		       &dimension,
		       this->m_P->GetX(),
		       this->m_P->GetB(),
		       this->m_P->GetI(),			   
		       tmp_,
		       y_);

#else
	  *outHasFailed_ = true;
#endif
	  
	};
	  
      };

    }; // namespace MKL

  }; // namespace Preconditioner
  
}; // namespace WLS

