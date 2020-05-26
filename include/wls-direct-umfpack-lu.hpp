#pragma once

#include "IInverseOperator.hpp"

#if __MK_WITH_UMFPACK__
#if __MK_ILP64__
#define UF_long long long
#define UF_long_max 9223372036854775801
#define UF_long_id "%Ld"
#define umfpack_symbolic 	umfpack_dl_symbolic
#define umfpack_solve    	umfpack_dl_solve
#define umfpack_numeric  	umfpack_dl_numeric
#define umfpack_symbolic_free   umfpack_dl_free_symbolic
#define umfpack_numeric_free    umfpack_dl_free_numeric
#else
#define umfpack_symbolic_free   umfpack_dl_free_symbolic
#define umfpack_numeric_free    umfpack_dl_free_numeric
#define umfpack_symbolic 	umfpack_dl_symbolic
#define umfpack_solve    	umfpack_dl_solve
#define umfpack_numeric  	umfpack_dl_numeric
#endif
#include "umfpack.h"
#endif

namespace WLS
{
  namespace Direct
  {
    namespace UMFPACK
    {
      
      /// <summary>
      /// Implementation of the linear solver with PARDISO.
      /// </summary>
      class Umfpack : public IInverseOperator
      {

      private:
	/// <summary>
	/// Enumerate the matrix types.
	/// </summary>
	struct MatrixType
	{
	public:
	  typedef enum _Value
	    {
	      /// <summary>
	      /// Real and structurally symmetric.
	      /// </summary>
	      RealAndStructurallySymmetric = 1,
	      
	      /// <summary>
	      /// Real and symmetric positive definite.
	      /// </summary>
	      RealAndSymmetricPositiveDefinite = 2,
	      
	      /// <summary>
	      /// Real and symmetric indefinite.
	      /// </summary>
	      RealAndSymmetricIndefinite = -2,
	      
	      /// <summary>
	      /// Real and not symmetric.
	      /// </summary>
	      ReadAndUnsymmetric = 11
	    } Value;
	};
	
	struct Error
	{
	public:
	  typedef enum _Value
	    {
	      /// <summary>
	      /// no error.
	      /// </summary>
	      OK = UMFPACK_OK,
	
	      /// <summary>
	      /// input inconsistency.
	      /// </summary>
	      invalid_matrix = UMFPACK_ERROR_invalid_matrix,

	      /// <summary>
	      /// not enough memory.
	      /// </summary>
	      invalid_Symbolic_object = UMFPACK_ERROR_invalid_Symbolic_object
	    } Value;
	};

    
	/// <summary>
	/// Enumerate options for the matrix reordering.
	/// </summary>
	struct MatrixReordering
	{
	public:
	  typedef enum _Value
	    {
	      /// <summary>
	      /// 
	      /// </summary>
	      MinimumDegree = 0,
	      /// <summary>
	      /// 
	      /// </summary>
	      NestedDissectionFromMetis = 2,
	      /// <summary>
	      /// 
	      /// </summary>
	      OpenMPVersionOfNestedDissection = 3
	    } Value;
	};
	

	struct Job
	{
	public:
	  typedef enum _Value
	    {
	      /// <summary>
	      /// Release the memory.
	      /// </summary>
	      ReleaseMemory = -1,
	
	      /// <summary>
	      /// Perform symbolic factorization.
	      /// </summary>
	      SymbolicFactorization = 11,
        
	      /// <summary>
	      /// Perform numerical factorization.
	      /// </summary>
	      NumericalFactorization = 22,
	
	      /// <summary>
	      /// Solve.
	      /// </summary>
	      Solve = 33,
	    } Value;
	};
	
  
	class Parameters
	{
	private:
	  /**
	     \brief Umfpack flag
	  */
	  Error::Value m_error;
	  /**
	     \brief Umfpack symbolic pointer
	  */
	  void * umfpack_Symbolic;
	  /**
	     \brief Umfpack numeric pointer
	  */
	  void * umfpack_Numeric;
	  /**
	     \brief Umfpack info
	  */
	  nsREAL umfpack_Info 		[UMFPACK_INFO];
	  /**
	     \brief Umfpack control
	  */
	  nsREAL umfpack_Control 	[UMFPACK_CONTROL];
	
	public:
	  inline void**GetHandle()
	  {
	    return m_pt; 
	  };
    
	  inline bool IsFortranIndexing() const
	  {
	    return this->m_iparm[34] == 0;      
	  };
    
	  inline void SetFortranIndexing(const bool value) 
	  {
	    this->m_iparm[34] = value ? 0 : 1;
	  };
    
	  inline bool IsTransposed()const
	  {
	    return this->m_iparm[11] == 2;       
	  };
    
	  inline void SetTransposed(const bool value)
	  {
	    this->m_iparm[11] = (value) ? 2 : 0;
	  };



	  /// <summary>
	  /// Get/Set the matrix reordering to apply.
	  /// </summary>
	  inline MatrixReordering::Value GetReordering()const
	  {
	    return (MatrixReordering::Value)m_iparm[1]; 
	  };

	  inline void SetReordering(const MatrixReordering::Value value)
	  {
	    this->m_iparm[1] = (int)value;
	  };
    
	  /// <summary>
	  /// The arrya of the integer parameters.
	  /// </summary>
	  inline WLS::integer_pt GetParamIntegers()
	  {
	    return this->m_iparm;
	  };

	  /// <summary>
	  /// Constructor.
	  /// </summary>
	  Parameters()
	  {
	    //
	    // Set the handle to zero.
	    //
	    //	  for (I i = 0; i < 64; i++)
	    //	    {
	    //	      m_pt[i] = NULL;
	    //	    }
	    //	  for (I i = 0; i < 64; i++)
	    //	    {
	    //	      m_iparm[i] = 0;
	    //	    }
	    //		
	    m_iparm[0] = 1; // No solver default
	    m_iparm[1] = 2; // Fill-in reordering from METIS 
	    // Numbers of processors, value of OMP_NUM_THREADS 
	    m_iparm[2] = 32;
	    m_iparm[3] = 31; // No iterative-direct algorithm 
	    m_iparm[4] = 0; // No user fill-in reducing permutation 
	    m_iparm[5] = 0; // Write solution into x 

	    //
	    // OUTPUT:NUmber of iterative refinement steps performed
	    //
	    m_iparm[6] = 0; // Not in use 


	    m_iparm[7] = 2; // Max numbers of iterative refinement steps 
	    m_iparm[8] = 0; // Not in use 
	    m_iparm[9] = 13; // Perturb the pivot elements with 1E-13 

	    //
	    // Scaling factors
	    //
	    m_iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS 

	    m_iparm[12] = 1;
	    // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy 

	    m_iparm[13] = 0; // Output: Number of perturbed pivots 
	    m_iparm[14] = 0; // Not in use 
	    m_iparm[15] = 0; // Not in use 
	    m_iparm[16] = 0; // Not in use 
	    m_iparm[17] = -1; // Output: Number of nonzeros in the factor LU 
	    m_iparm[18] = -1; // Output: Mflops for LU factorization 
	    m_iparm[19] = 0; // Output: Numbers of CG Iterations 

	    //
	    // We choose C-indexing by default
	    //
	    this->SetFortranIndexing(false);

	    //
	    // We solve Ax=b by default
	    //
	    this->SetTransposed(false);

	    //
	    // Reordering
	    //
	    this->SetReordering(MatrixReordering::MinimumDegree);

	  };


    
	};
	
	
      private:
	/// <summary>
	/// Type of the matrix.
	/// </summary>
	static const WLS::integer_t s_mtype = (int) MatrixType::ReadAndUnsymmetric;
	
	/// <summary>
	/// Maximal number of factors in memory. Generally used value is 1.
	/// </summary>
	static const WLS::integer_t s_maxfct = 1;

	/// <summary>
	/// The number of matrix (from 1 to <see cref="s_maxfct"/>) to solve.
	/// </summary>
	static const WLS::integer_t s_mnum = 1;
	
	/// <summary>
	/// Message level information.
	/// </summary>
	/// <remarks>
	/// 0: PARDISO generates no output
	/// 1: PARDISO prints statistical information
	/// </remarks>
	static const WLS::integer_t s_msglvl = 0;
	
	/// <summary>
	/// We fix the number of right-hand sides.
	/// </summary>
	static const WLS::integer_t s_nrhs = 1;

	/// <summary>
	/// Array of pointers.
	/// </summary>
	Parameters m_parameters;

	/// <summary>
	/// The error variable.
	/// </summary>
	Error::Value m_error;
	
	/// <summary>
	/// The matrix.
	/// </summary>

	WLS::integer_t  m_size{0};
	WLS::integer_pt m_begin   {NULL};
	WLS::integer_pt m_indices {NULL};
	double*__restrict__ m_values  {NULL};
	
	bool m_symbolicFactorizationDone {false};
	
      private:
	struct Error
	{
	public:
	  typedef enum _Value
	    {
	      OK = UMFPACK_OK,
	      warning_singular_matrix = UMFPACK_WARNING_singular_matrix,
	      warning_determinant_underflow = UMFPACK_WARNING_determinant_underflow,
	      warning_determinant_overflow = UMFPACK_WARNING_determinant_overflow
	      out_of_memory = UMFPACK_ERROR_out_of_memory,
	      invalid_Numeric_object = UMFPACK_ERROR_invalid_Numeric_object,
	      invalid_Symbolic_object = UMFPACK_ERROR_invalid_Symbolic_object,
	      argument_missing = UMFPACK_ERROR_argument_missing,	      
	      n_nonpositive =  UMFPACK_ERROR_n_nonpositive,
	      invalid_matrix = UMFPACK_ERROR_invalid_matrix,
	      different_pattern = UMFPACK_ERROR_different_pattern,
	      invalid_system = UMFPACK_ERROR_invalid_system,
	      invalid_permutation = UMFPACK_ERROR_invalid_permutation,	      
	      file_IO = UMFPACK_ERROR_file_IO,	      
	      ordering_failed = UMFPACK_ERROR_ordering_failed,
	      internal_error = UMFPACK_ERROR_internal_error
	    } Value;
	};
	
	static std::string ErrorMessage(const Error::Value error)
	{
	  switch(error)
	    {
	    case Error::OK:
	      {
		return "UMFPACK was successful.";
	      }
	    case Error::warning_singular_matrix:
	      {
		return "Matrix is singular. There are exact zeros on the diagonal of U.";
	      }
	      
	    case Error::warning_determinant_underflow:
	      {
		return "The determinant is nonzero, but smaller in magnitude than the smallest positive floating-point number.";
	      }
	    case Error::warning_determinant_overflow:
	      {
		return "The determinant is larger in magnitude than the largest positive floating-point number (IEEE Inf).";
	      }
	    case Error::out_of_memory:
	      {
		return "Not enough memory. The ANSI C malloc or realloc routine failed.";
	      }
	    case Error::invalid_Numeric_object:
	      {
		return "Routines that take a Numeric object as input (or load it from a file) check this object and return this error code if it is invalid."
		  " This can be caused by a memory leak or overrun in your program, which can overwrite part of the Numeric object."
		  " It can also be caused by passing a Symbolic object by mistake, or some other pointer."
		  " If you try to factorize a matrix using one version of UMFPACK and then use the factors in another version, this error code will trigger as well. "
		  "You cannot factor your matrix using version 4.0 and then solve with version 4.1, for example. 3 . You cannot use different precisions of the same version (real and complex, for example). It is possible for the Numeric object to be corrupted by your program in subtle ways that are not detectable by this quick check. In this case, you may see an UMFPACK ERROR different pattern error code, or even an UMFPACK ERROR internal error.";
	      }
	    case Error::invalid_Symbolic_object:
	      {
		return "Routines that take a Symbolic object as input (or load it from a file) check this object and return this error code if it is invalid. The causes of this error are analogous to the UMFPACK ERROR invalid Numeric object error described above.";
	      }
	    case Error::argument_missing:
	      {
		return "Some arguments of some are optional (you can pass a NULL pointer instead of an array). This error code occurs if you pass a NULL pointer when that argument is required to be present.";
	      }
	      
	    case  Error::n_nonpositive:
	      {
		
		return "The number of rows or columns of the matrix must be greater than zero.";
	      }																					case Error::invalid_matrix:
	      {
		return "The matrix is invalid. For the column-oriented input, this error code will occur if the contents of Ap and/or Ai are invalid. Ap is an integer array of size n col+1. On input, it holds the “pointers” for the column form of the sparse matrix A. Column j of the matrix A is held in Ai [(Ap [j]) . . . (Ap [j+1]-1)]. The first entry, Ap [0], must be zero, and Ap [j] ≤ Ap [j+1] must hold for all j in the range 0 to n col-1. The value nz = Ap [n col] is thus the total number of entries in the pattern of the matrix A. nz must be greater than or equal to zero. The nonzero pattern (row indices) for column j is stored in Ai [(Ap [j]) . . . (Ap [j+1]-1)]. The row indices in a given column j must be in ascending order, and no duplicate row indices may be present. Row indices must be in the range 0 to n row-1 (the matrix is 0-based). Some routines take a triplet-form input, with arguments nz, Ti, and Tj. This error code is returned if nz is less than zero, if any row index in Ti is outside the range 0 to n col-1, or if any column index in Tj is outside the range 0 to n row-1.";
	      }
	    case Error::different_pattern:
	      {
		return "The most common cause of this error is that the pattern of the matrix has changed between the symbolic and numeric factorization. It can also occur if the Numeric or Symbolic object has been subtly corrupted by your program.";
	      }
	    case Error::invalid_system:
	      {
		return "The sys argument provided to one of the solve routines is invalid.";
	      }
	    case Error::invalid_permutation:
	      {
		return "The permutation vector provided as input is invalid.";
	      }
	      
	    case Error::file_IO:
	      {
		return "This error code is returned by the routines that save and load the Numeric or Symbolic objects to/from a file, if a file I/O error has occurred. The file may not exist or may not be readable, you may be trying to create a file that you don’t have permission to create, or you may be out of disk space. The file you are trying to read might be the wrong one, and an earlier end-of-file condition would then result in this error.";
	      }
	    case Error::ordering_failed:
	      {
		return "The ordering method failed.";
	      }
	    case Error::internal_error:
	      {
		return "An internal error has occurred, of unknown cause. This is either a bug in UMFPACK, or the result of a memory overrun from your program. Try modifying the file AMD/Include/amd internal.h and adding the statement #undef NDEBUG, to enable the debugging mode. Recompile UMFPACK and rerun your program. A failed assertion might occur which can give you a better indication as to what is going wrong. Be aware that UMFPACK will be extraordinarily slow when running in debug mode. If all else fails, contact the developer (DrTimothyAldenDavis@gmail.com) with as many details as possible.";
	      }
	    default:
	      {
		return "Unknown value."
	      }
	    }			       
	  
	};

      
      static void CallUmfpack(Parameters *parameters,
			      WLS::integer_t maxfct,
			      WLS::integer_t mnum,
			      WLS::integer_t mtype,
			      const Job::Value job,
			      WLS::integer_t n,
			      double* __restrict__ a,
			      WLS::integer_pt ia,
			      WLS::integer_pt ja,
			      WLS::integer_pt perm,
			      WLS::integer_t nrhs,
			      WLS::integer_t msglvl,
			      double* __restrict__ b,
			      double* __restrict__ x,
			      Error::Value* outError)
      {
	WLS::integer_t err;
	WLS::integer_t phase = (WLS::integer_t) job;
	
	pardiso(parameters->GetHandle(),
		&maxfct,
		&mnum,
		&mtype,
		&phase,
		&n,
		a,
		ia,
		ja,
		perm,
		&nrhs,
		parameters->GetParamIntegers(),
		&msglvl,
		b,
		x,
		&err);
	
	*outError = (Error::Value) err;
      };
	
      public:	
	/// <summary>
	/// Constructor.
	/// </summary>
	/// <param name="A">The sparse matrix  <see cref="SparseMatrixCSR32"/>.</param>
	Umfpack(const WLS::integer_t 	size_,
		WLS::integer_pt 	begin_,
		WLS::integer_pt 	indices_,
		double* __restrict__ 	values_,
		bool fortranIndexing = true)
	  : m_size(size_),
	    m_begin(begin_),
	    m_indices(indices_),
	    m_values(values_)
	{ 
	  this->m_error = Error::NoError;	  
	  this->m_parameters.SetFortranIndexing(fortranIndexing);
	  this->m_symbolicFactorizationDone = false;
	  umfpack_dl_defaults (umfpack_Control);	
	};
	
	void Apply(const char * transpose,
		   double* __restrict__ y,
		   const double* __restrict__  x,
		   WLS::integer_t tmpSize,
		   double* __restrict__ tmp,
		   bool* outHasFailed)
	{

#if __MK_WITH_UMFPACK__
	  umfpack_solve(UMFPACK_A,
			this->m_begin,
			this->m_indices,
			this->m_values,
			y,
			x,
			this->umfpack_Numeric,
			this->umfpack_Control,
			this->umfpack_Info);
#else
	  std::cerr << "mkFACTORIZATION_lusol:umfpack library not available" << std::endl;
#endif
	  *outHasFailed = m_error != Error::NoError;
	  
	  //
	  // We do not care about the temporary vector tmp.
	  //
	  //	  m_parameters.SetTransposed(transpose[0] == 'T' || transpose[0] == 't');
      }

	
	/// <see cref="IInverseOperator.Compute"/>
	void Compute(bool* outHasFailed)
	{
#if __MK_WITH_UMFPACK__  
	    if (!this->m_symbolicFactorizationDone)
	      {	    
		if (!this->umfpack_Symbolic)
		  {
		    this->m_error = (Error::Value)umfpack_symbolic(N,
								   N, 
								   m_begin,
								   m_indices,
								   m_values,
								   &this->umfpack_Symbolic,
								   this->umfpack_Control,
								   this->umfpack_Info);
		    
		    if (*outHasFailed = (this->m_error != Error::OK ) )
		      {
			return;
		      }
		  }
		
		this->m_symbolicFactorizationDone = true;
	      }


	    if (this->umfpack_Numeric)
	      {
		umfpack_numeric_free(&this->umfpack_Numeric);
		this->umfpack_Numeric = nullptr;
	      }

	    this->m_error = (Error::Value)umfpack_numeric(m_begin,
							  m_indices,
							  m_value,
							  this->umfpack_Symbolic, 
							  &this->umfpack_Numeric, 
							  this->umfpack_Control,
							  this->umfpack_Info);
	    if (*outHasFailed = (this->m_error != Error::OK ) )
	      {
		return;
	      }
	    
	};

	virtual WLS::integer_t GetSizeOfTemporaryVector() const
	{
	  return 0;
	};
	

	virtual std::string GetErrorMessage() const
	{
	  return ErrorMessage(m_error); 
	};
	

	~Umfpack()
	{

	  if (this->umfpack_Symbolic)
	    {
	      this->umfpack_symbolic_free(&this->umfpack_Symbolic);
	      this->umfpack_Symbolic=NULL;
	    }
	  
	  if (this->umfpack_Numeric)
	    {
	      this->umfpack_numeric_free(&this->umfpack_Numeric);
	      this->umfpack_Numeric=NULL;
	    }
 
	};
	
      };

    }; // namespace MKL
    
  }; // namespace Direct
  
}; // namespace WLS
