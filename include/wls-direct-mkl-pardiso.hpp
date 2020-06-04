#pragma once
#include "wls.hpp"

namespace WLS
{
  namespace Direct
  {
    namespace MKL
    {
      
      /// <summary>
      /// Implementation of the linear solver with PARDISO.
      /// </summary>
      class Pardiso : public inverse_operator
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
	      NoError = 0,
	
	      /// <summary>
	      /// input inconsistency.
	      /// </summary>
	      InconsistentInput = -1,

	      /// <summary>
	      /// not enough memory.
	      /// </summary>
	      NotEnoughMemory = -2,
               
	      /// <summary>
	      /// reordering problem.
	      /// </summary>
	      ReorderingProblem = -3,

	      /// <summary>
	      /// zero pivot, numerical factorization or iterative refinement problem.
	      /// </summary>
	      ZeroPivotNumericalFactorizationOrRefinementProblem = -4,

	      /// <summary>
	      /// unclassified (internal) error.
	      /// </summary>
	      Unclassified = -5,

	      /// <summary>
	      /// pre-ordering failed (matrix types 11, 13 only).
	      /// </summary>
	      PreorderingFailed = -6,

	      /// <summary>
	      /// diagonal matrix is singular.
	      /// </summary>
	      SingularDiagonalMatrix = -7,

	      /// <summary>
	      /// 32-bit integer overflow problem.
	      /// </summary>
	      Int32Overflow = -8,

	      /// <summary>
	      /// not enough memory for OOC.
	      /// </summary>
	      NotEnoughMemoryForOOC = -9,
               
	      /// <summary>
	      /// problems with opening OOC temporary file.
	      /// </summary>
	      ProblemsWithOpeningOOCTemporaryFiles = -10,

	      /// <summary>
	      /// read/write problems with the OOC data file.
	      /// </summary>
	      ReadWriteProblemsWithTheOOCDataFile = -11
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
	  wls_int_t	m_iparm[64] {0};
	  void * 	m_pt[64] {NULL};
	
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
	  inline wls_int_p GetParamIntegers()
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
	//!
	//! @brief Type of the matrix.
	//!
	static const wls_int_t s_mtype = (int) MatrixType::ReadAndUnsymmetric;
	
	//!
	//! @brief Maximal number of factors in memory. Generally used value is 1.
	//!
	static const wls_int_t s_maxfct = 1;

	//!
	//! @brief The number of matrix (from 1 to <see cref="s_maxfct"/>) to solve.
	//!
	static const wls_int_t s_mnum = 1;
	
	//!
	//! @brief Message level information.
	//!
	/// <remarks>
	/// 0: PARDISO generates no output
	/// 1: PARDISO prints statistical information
	/// </remarks>
	static const wls_int_t s_msglvl = 0;
	
	//!
	//! @brief We fix the number of right-hand sides.
	//!
	static const wls_int_t s_nrhs = 1;

	//!
	//! @brief Array of pointers.
	//!
	Parameters m_parameters;

	//!
	//! @brief The error variable.
	//!
	Error::Value m_error;
	
	//! @brief The matrix.

	wls_int_t  	m_size		{0};
	wls_int_p 	m_begin   	{NULL};
	wls_int_p 	m_indices 	{NULL};
	double*__restrict__ m_values  		{NULL};
	
	bool m_symbolicFactorizationDone 	{false};
	
      private:

	static std::string ErrorMessage(const Error::Value error)
	{	
	  switch (error)
	    {
	    case Error::InconsistentInput:
	      {
		return "input inconsistent";
	      }
	    case Error::NotEnoughMemory:
	      {
		return "not enough memory";
	      }
	    case Error::ReorderingProblem:
	      {
		return "reordering problem";
	      }
	    case Error::ZeroPivotNumericalFactorizationOrRefinementProblem:
	      {
		return "zero pivot, numerical factorization or iterative refinement problem";
	      }
	    case Error::Unclassified:
	      {
		return "unclassified (internal) error";
	      }
	    case Error::PreorderingFailed:
	      {
		return "pre-ordering failed (matrix types 11, 13 only)";
	      }
	    case Error::SingularDiagonalMatrix:
	      {
		return "diagonal matrix is singular";
	      }
	    case Error::Int32Overflow:
	      {
		return "32-bit integer overflow problem";
	      }
	    case Error::NotEnoughMemoryForOOC:
	      {
		return "not enough memory for OOC";
	      }
	    case Error::ProblemsWithOpeningOOCTemporaryFiles:
	      {
		return "problems with opening OOC temporary file";
	      }
	    case Error::ReadWriteProblemsWithTheOOCDataFile:
	      {
		return "read/write problems with the OOC data file";
	      }
	    
	    default:
	      {
#if 0
		Debug.Assert(Error.NoError == error, "Unexpected value of the enumeration Error.");
#endif
		return "no error";
	      }
	    }
	};

      
      static void CallPardiso(Parameters *parameters,
			      wls_int_t maxfct,
			      wls_int_t mnum,
			      wls_int_t mtype,
			      const Job::Value job,
			      wls_int_t n,
			      double* __restrict__ a,
			      wls_int_p ia,
			      wls_int_p ja,
			      wls_int_p perm,
			      wls_int_t nrhs,
			      wls_int_t msglvl,
			      double* __restrict__ b,
			      double* __restrict__ x,
			      Error::Value* outError)
      {
#ifdef WLS_WITH_MKL
	wls_int_t err;
	wls_int_t phase = (wls_int_t) job;
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
#endif
      };
	
      public:	
	Pardiso(const wls_int_t 	size_,
		wls_int_p 	begin_,
		wls_int_p 	indices_,
		double* __restrict__ 	values_,
		bool fortranIndexing = true)
	  : m_size(size_),
	    m_begin(begin_),
	    m_indices(indices_),
	    m_values(values_)
	{ 
	  this->m_error = Error::NoError;	  
	  //
	  // Ensure we have C-indexing
	  //
	  this->m_parameters.SetFortranIndexing(fortranIndexing);
	  this->m_symbolicFactorizationDone = false;
	};
	


	void apply(const char * transpose,
		   double* __restrict__ y,
		   const double* __restrict__  x,
		   wls_int_t tmpSize,
		   double* __restrict__ tmp,
		   bool* outHasFailed)
	{
	  
         //
         // We do not care about the temporary vector tmp.
         //

	  m_parameters.SetTransposed(transpose[0] == 'T' || transpose[0] == 't');
	  CallPardiso(&m_parameters,
		      s_maxfct,
		      s_mnum,
		      s_mtype,
		      Job::Solve,		      
		      m_size,
		      m_values,
		      m_begin,
		      m_indices,			     		      
		      NULL, // perm
		      s_nrhs,
		      s_msglvl,
		      (double* __restrict__)x, // b
		      y, // x
		      &this->m_error);
	  *outHasFailed = m_error != Error::NoError;
	  
      }


	void compute(bool* outHasFailed)
	{
	  if (!this->m_symbolicFactorizationDone)
	    {
	      
	      this->m_symbolicFactorizationDone = true;
	      //
	      // Compute the symbolic factorization.
	      // 
	      CallPardiso(&m_parameters,
			  s_maxfct,
			  s_mnum,
			  s_mtype,
			  Job::SymbolicFactorization,
			  m_size,
			  m_values,
			  m_begin,
			  m_indices,			     
			  NULL,//perm
			  s_nrhs,
			  s_msglvl,
			  NULL, // b
			  NULL, // x
			  &m_error);

	      if (m_error != Error::NoError)
		{
		  *outHasFailed = true;
		  return;
		}	  

	    }
	  
	  CallPardiso(&m_parameters,
		      s_maxfct,
		      s_mnum,
		      s_mtype,
		      Job::NumericalFactorization,
		      m_size,
		      m_values,
		      m_begin,
		      m_indices,			     		      		      
		      NULL, // perm
		      s_nrhs, 
		      s_msglvl,
		      NULL, // b
		      NULL, // x
		      &this->m_error);

         *outHasFailed = m_error != Error::NoError;
	};

	virtual wls_int_t get_buffer_size() const
	{
	  return 0;
	};
	
	virtual std::string get_error_message() const
	{
	  return ErrorMessage(m_error); 
	};
	
	virtual ~Pardiso()
	{
	  
	  CallPardiso(&this->m_parameters,
		      s_maxfct,
		      s_mnum,
		      s_mtype,
		      Job::ReleaseMemory,
		      m_size,
		      NULL,
                      m_begin,
		      m_indices,		      
                      NULL,
                      s_nrhs,
                      s_msglvl,
                      NULL,
                      NULL,
		      &this->m_error);

	};
	
      };

    }; // namespace MKL
    
  }; // namespace Direct
  
}; // namespace WLS
