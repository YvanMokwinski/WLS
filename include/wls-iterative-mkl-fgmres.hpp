#pragma once

#include "ILinearOperator.hpp"
#include "Iterative/MKL/Parameters.hpp"
#include "WLA/include/BlasMKL.hpp"
#include "WLS_MKL.hpp"

namespace WLS
{
  namespace Iterative
  {
    namespace MKL
    {
      //!
      //! @brief Implementation of the Flexible General Minimal Residual with the Math Kernel Library.
      //!
      class Fgmres : public WLS::IInverseOperator
      {
 
      private:

	struct Status
	{
	public:
	  //!
	  //! @brief Enumerate the values of the parameter request.
	  //!
	  typedef enum _Value
	    {
	      //!
	      //! @brief Indicates that the task completed normally and the solution is found and stored in the vector x. 
	      //! This occurs only if the stopping tests are fully automatic. 
	      //! For the user defined stopping tests, see the description of the RCI_request= 2 or 4.
	      //!
	      Completed = 0,
	  
	      //!
	      //! @brief Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met. 
	      //! This situation occurs only if you request both tests.
	      //!
	      MaximumNumberOfIterationsReached = -1,
	  
	      //!
	      //! @brief Indicates that the routine was interrupted because of an attempt to divide by zero. 
	      //!  Usually this happens if the matrix is degenerate or almost degenerate. 
	      //! However, it may happen if the parameter dpar is altered, or if the method is not stopped when the solution is found.
	      //!
	      AttemptDivideByZero = -10,
	  
	      //!
	      //! @brief Indicates that the routine was interrupted because it entered an infinite cycle. 
	      //! Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine, 
	      //! or the dfgmres_check routine was not called.
	      //!
	      InfiniteCycle = -11,

	      //!
	      //! @brief Indicates that the routine was interrupted because errors were found in the method parameters. 
	      //! Usually this happens if the parameters ipar and dpar were altered by mistake outside the routine.
	      //!
	      InconsitentParameters = -12,

	      //!
	      //! @brief Indicates that you must multiply the matrix by rwork[ipar[21] - 1:ipar[21] + n - 2], 
	      //! put the result in the rwork[ipar[22] - 1:ipar[22] + n - 2], and return control back to the routine dfgmres.
	      //!
	      ApplyMatrixVectorProduct = 1,

	      //!
	      //! @brief Indicates that you must perform the stopping tests. If they fail, return control to the dfgmres routine. Otherwise, the FGMRES solution is found, and you can run the fgmres_get routine to update the computed solution in the vector x.
	      //!
	      ApplyStoppingTests = 2,

	      //!
	      //! @brief Indicates that you must apply the inverse preconditioner to tmp[ipar[21] - 1:ipar[21] + n - 2], put the result in the tmp[ipar[22] - 1:ipar[22] + n - 2], and return control back to the routine dfgmres.
	      //!
	      ApplyPreconditioner = 3,

	      //!
	      //! @brief Indicates that you must check the norm of the currently generated vector. 
	      //! If it is not zero within the computational/rounding errors, return control to the dfgmres routine. 
	      //! Otherwise, the FGMRES solution is found, and you can run the dfgmres_get routine to update the computed solution in the vector x.
	      //!
	      CheckTheNormOfTheGeneratedVector = 4
	    } Value;


	  //!
	  //! @brief Get the error message related to the status.
	  //!
	  //! @param request The parameter 'request' from the routine.
	  //!
	  static std::string GetErrorMessage(const Status::Value request)
	  {
	    switch (request)
	      {
	      case MaximumNumberOfIterationsReached:
		{
		  return
		    "Indicates that the routine was interrupted because the maximum number of iterations was reached, but the relative stopping criterion was not met. This situation occurs only if you request both tests.";
		}
	     
	      case AttemptDivideByZero:
		{
		  return
		    "Indicates that the routine was interrupted because of an attempt to divide by zero. Usually this happens if the matrix is degenerate or almost degenerate. However, it may happen if the parameter dpar is altered, or if the method is not stopped when the solution is found.";
		}
	    
	      case InfiniteCycle:
		{
		  return "Indicates that the routine was interrupted because it entered an infinite cycle. Usually this happens because the values ipar[7], ipar[8], ipar[9] were altered outside of the routine, or the dfgmres_check routine was not called.";
		}
	    
	      case InconsitentParameters:
		{
		  return "Indicates that the routine was interrupted because errors were found in the method parameters. Usually this happens if the parameters ipar and dpar were altered by mistake outside the routine.";
		}
	    
	      case Completed:
		{
		  return "Completed";
		}
	      case ApplyMatrixVectorProduct:
		{
		  return "ApplyMatrixVectorProduct";
		}
	      case ApplyStoppingTests:
		{
		  return "ApplyStoppingTests";
		}
	      case ApplyPreconditioner:
		{
		  return "ApplyPreconditioner";
		}
	      case CheckTheNormOfTheGeneratedVector:
		{
		  return "ApplyPreconditioner";
		}
	      default:
		{
		  int irequest = (int) request;
		  if (irequest < 0)
		    {
		      switch (irequest)
			{
			case -10000:
			  {
			    return "Indicates that the routine failed to complete the task.";
			  }
			default:
			  {
			    return "Unknown error";
			  }
			}
		    }
		  else
		    {
		      return "Unknown request value";
		    }
		}
	      }
	  };
       
	};

      private:
          
     
	//!
	//! @brief Parameters of the MKL FGMRES Solver.
	//!
	class ParametersFgmres : public Parameters
	{
            
	public:
      
	  ParametersFgmres()
	  {
	  };

	public:
	  virtual ~ParametersFgmres()
	  {
	  };

	public:
       
	  //!
	  //! @brief Set the relative tolerance.
	  //!
	  void SetRelativeTolerance(const double value)
	  {
	    this->m_dpar[0] = value; 
	  };

	  //!
	  //! @brief Make the residual stopping test automatic.
	  //!
	  void SetAutomaticResidualStoppingTest(const bool value)
	  {
	    this->m_ipar[8] = value ? 1 : 0;
	  };
	
	  //!
	  //! @brief Set the maximum number of iterations.
	  //!
	  void SetMaximumNumberOfIterations(const wls_int_t  value)
	  {
	    this->m_ipar[4] = value;
	  }
	
	  //!
	  //! @brief Set the maximum number of non-restarted iterations.
	  //!
	  void SetMaximumNumberOfNonRestartedIterations(const wls_int_t  value)
	  {
	    this->m_ipar[14] = value;
	  };
        
	
	  //!
	  //! @brief Get the norm of the orthogonal vector.
	  //!
	  double GetNormOfOrthogonalVector()const
	  {
	    return this->m_dpar[6]; 
	  };

	  //!
	  //! @brief Get/Set the threshold of the norm of the orthogonal vector.
	  //!
	  double GetZeroNormThresholdOfOrthogonalVector() const
	  {
	    return this->m_dpar[7]; 
	  };

      
	  void SetZeroNormThresholdOfOrthogonalVector(const double value) 
	  {
	    this->m_dpar[7] = value;
	  };
	
	  //!
	  //! @brief Get the index of the source temporary vector.
	  //!
	  wls_int_t  GetIndexOfSourceTemporaryVector()const
	  {
	    return this->m_ipar[21] - 1;
	  };
	
	  //!
	  //! @brief Get the index of output temporary vector.
	  //!
	  wls_int_t  GetIndexOfOutputTemporaryVector()const
	  {
	    return this->m_ipar[22] - 1;
	  };
	
	  //!
	  //! @brief Make the user defined stopping test automatic.
	  //!
	  void SetAutomaticUserDefinedStoppingTest(const bool value)
	  {
	    this->m_ipar[9] = value ? 1 : 0; 
	  };
	
	  //!
	  //! @brief Make the zero norm of the orthogonal vector automatic.
	  //!
	  bool GetAutomaticTestForZeroNormOfOrthogonalVector()
	  {
	    return this->m_ipar[11] != 0; 	  
	  };

	  void SetAutomaticTestForZeroNormOfOrthogonalVector(const bool value)
	  {
	    this->m_ipar[11] = value ? 1 : 0;
	  };
	
	  //!
	  //! @brief Indicates if the method is preconditioned.
	  //!
	  bool IsPreconditionedMethod()const
	  {
	    return this->m_ipar[10] == 1; 
	  };
	 
	  void SetPreconditionedMethod(const bool value)
	  {	  
	    this->m_ipar[10] = value ? 1 : 0; 
	  };
       
	  //!
	  //! @brief Print information.
	  //!
	  //! @param textWriter The text writer.
	  void PrintInformation()
	  {
	    std::cout << "Some info about the current run of RCI FGMRES method:"<< std::endl;
	    if (m_ipar[7] != 0)
	      {
		std::cout << "   As ipar[7]="<< m_ipar[7]<< ", the automatic test for the maximal number of iterations will be performed"<< std::endl;
	      }
	    else
	      {
		std::cout << "   As ipar[7]="<< m_ipar[7]<< ", the automatic test for the maximal number of iterations will be skipped"<< std::endl;
	      }
	    std::cout << "+++"<< std::endl;
	    if (0 != m_ipar[8])
	      {
		std::cout << "   As ipar[8]="<< m_ipar[8]<< ", the automatic residual test will be performed" << std::endl;
	      }
	    else
	      {
		std::cout << "   As ipar[8]="<< m_ipar[8]<< ", the automatic residual test will be skipped" << std::endl;
	      }
	    std::cout << "+++"<< std::endl;
	    if (0 != m_ipar[9])
	      {
		std::cout << "   As ipar[9]="<< m_ipar[9]<< ", the user-defined stopping test will be requested via RCI_request=2"<< std::endl;
	      }
	    else
	      {
		std::cout << "   As ipar[9]="<< m_ipar[9]<< ", the user-defined stopping test will not be requested, thus RCI_request will not take the value 2"<< std::endl;
	      }
	    std::cout << "+++"<< std::endl;
	    if (0 != m_ipar[10])
	      {
		std::cout << "   As ipar[10]="<< m_ipar[10]<< ", the Preconditioned FGMRES iterations will be performed, thus the preconditioner action will be requested via RCI_request=3"<< std::endl;
	      }
	    else
	      {
		std::cout << "   As ipar[10]="<< m_ipar[10]<< ", the Preconditioned FGMRES iterations will not be performed, thus RCI_request will not take the value 3"<< std::endl;
	      }
	    std::cout << "+++"<< std::endl;
	    if (0 != m_ipar[11])
	      {
		std::cout << "   As ipar[11]="<< m_ipar[11]<< ", the automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors will be performed, thus RCI_request will not take the value 4"<< std::endl;
	      }
	    else
	      {
		std::cout << "   As ipar[11]="<< m_ipar[11]<< ", the automatic test for the norm of the next generated vector is not equal to zero up to rounding and computational errors will be skipped, thus the user-defined test will be requested via RCI_request=4"<< std::endl;
	      }
	    std::cout << "+++"<< std::endl;
	  };

	  //!
	  //! @brief The required memory.
	  //!
	  //! @param maximumNumberOfNonRestartedIterations The maximum number of non-restarted iterations.
	  //! @param n The dimension.
	  /// @return The required memory.
	  static wls_int_t  GetRequiredMemory(const int maximumNumberOfNonRestartedIterations, const wls_int_t  n)
	  {
	    return ((2 * maximumNumberOfNonRestartedIterations + 1) * n +
		    maximumNumberOfNonRestartedIterations * (maximumNumberOfNonRestartedIterations + 9) / 2 + 1);
	  };

       
	};

     
      private:
     
	//!
	//! @brief Initializes the solver.
	//!
	//! @param n INPUT:Sets the size of the problem.
	//! @param request OUTPUT:Gives information about the result of the routine.
	//! @param parameters OUTPUT: The parameters. 
	//! @param rwork OUTPUT: Array of size ((2*parameters[14] + 1)*n + parameters[14]*(parameters[14 - 1] + 9)/2 + 1).
	static void Init(const wls_int_t  		n,
			 Status::Value* 	refStatus,
			 ParametersFgmres& 	parameters,
			 double*__restrict__ 		rwork)
	{
#ifdef WLS_WITH_MKL
	  wls_int_t  ln = n;
	  wls_int_t  status = (wls_int_t)*refStatus;

	  dfgmres_init(&ln,
		       NULL,
		       NULL,
		       &status,
		       parameters.GetParamIntegers(),
		       parameters.GetParamReals(),
		       rwork);
	  *refStatus = (Status::Value)status;
#endif
	};


	//!
	//! @brief Checks consistency and correctness of the user defined data.
	//!
	//! @param parameters OUTPUT: The parameters. 
	//! @param n INPUT:Sets the size of the problem.
	//! @param rwork OUTPUT: Array of size ((2*parameters [14] + 1)*n  + parameters [14]*(parameters [14] + 9)/2 + 1).
	//! @param request OUTPUT: Gives information about result of the routine.
	static void Check(ParametersFgmres&parameters,
			  const wls_int_t  n,
			  double*__restrict__ rwork,
			  Status::Value* refStatus)
	{
#ifdef WLS_WITH_MKL
	  wls_int_t  ln = n;
	  wls_int_t  status = (wls_int_t) *refStatus;
	  dfgmres_check(&ln,
			NULL,
			NULL,
			&status,
			parameters.GetParamIntegers(),
			parameters.GetParamReals(),
			rwork);
	  *refStatus = (Status::Value) status;
#endif
	};


	//!
	//! @brief Retrieves the number of the current iteration and updates the solution.
	//!
	//! @param parameters OUTPUT: The parameters. 
	//! @param n INPUT:Sets the size of the problem.
	//! @param sol OUTPUT: Array of size n . Contains the initial approximation to the solution vector. Normally it is equal to 0 or to rhs .
	//! @param rhs OUTPUT: Array of size n . Contains the right-hand side vector.
	//! @param request OUTPUT: Gives information about result of the routine.
	//! @param rwork INPUT:Array of size ((2*parameters [14] + 1)*n  + parameters [14]*(parameters [14] + 9)/2 + 1).
	//! @param itercount OUTPUT: Contains the value of the current iteration number.
	static void Get(ParametersFgmres&parameters,
			const wls_int_t  n,
			double*__restrict__ sol,
			double*__restrict__ rhs,
			double*__restrict__ rwork,
			Status::Value* refStatus,
			wls_int_t& itercount)
	{
#ifdef WLS_WITH_MKL	  
	  wls_int_t  ln = n;
	  wls_int_t  status = (wls_int_t)*refStatus;

	  dfgmres_get(&ln,
		      sol,
		      rhs,
		      &status,
		      parameters.GetParamIntegers(),
		      parameters.GetParamReals(),
		      rwork,
		      &itercount);
	  
	  *refStatus = (Status::Value) status;	
#endif
	};


	//!
	//! @brief Makes the FGMRES iterations.
	//!
	//! @param parameters OUTPUT: The parameters. 
	//! @param n INPUT:Sets the size of the problem.
	//! @param sol INPUT: Array of size n . Contains the initial approximation to the solution vector. Normally it is equal to 0 or to b.
	//! @param rhs INPUT: Array of size n . Contains the right-hand side vector.
	//! @param rwork INPUT/OUTPUT:The working array.
	//! @param request OUTPUT: Gives information about result of the routine.
	static inline void Run(ParametersFgmres& parameters,
			       const wls_int_t  n,
			       double*__restrict__ sol,
			       const double*__restrict__ rhs,
			       double*__restrict__ rwork,
			       Status::Value* refStatus)
	{	
#ifdef WLS_WITH_MKL
	  wls_int_t  ln = n;
	  wls_int_t  status = (wls_int_t) *refStatus;
	  dfgmres(&ln,
		  sol,
		  (double*__restrict__)rhs,
		  &status,
		  parameters.GetParamIntegers(),
		  parameters.GetParamReals(),
		  rwork);
	  *refStatus = (Status::Value) status;
#endif
	};
     
     
      private:
	//!
	//! @brief Operator for the matrix vector product.
	//!
	ILinearOperator* m_matrixVectorProductOperator;
     
	//!
	//! @brief The size of the linear system.
	//!
	wls_int_t  m_n;
     
	//!
	//! @brief The parameters.
	//!
	ParametersFgmres m_parameters;
     
	//!
	//! @brief The required memory for MKL.
	//!
	double*__restrict__ m_rwork;
     
	//!
	//! @brief The preconditioner.
	//!
	IInverseOperator* m_preconditioner;
     
	//!
	//! @brief Flag to indicate if any action of the preconditioner has failed.
	//!
	bool m_preconditionerHasFailed;
     
	//!
	//! @brief Status.
	//!
	Status::Value m_status;
     
	//!
	//! @brief The error message related to the status.
	//!
	std::string m_errorMessage;

     
	static const wls_int_t  s_maximumNumberOfNonRestartedIterations = 30;

	//!
	//! @brief Get the required memory to run the MKL Conjugate Gradient method. 
	//!
	//! @param size The size of the linear system to solve.
	//! @param maximumNumberOfNonRestartedIterations The maximum number of non-restarted iterations.
	/// @return The size of the required array of double.
	static wls_int_t  GetRequiredMemory(wls_int_t  size,
						  wls_int_t  maximumNumberOfNonRestartedIterations = s_maximumNumberOfNonRestartedIterations)
	{
	  return ParametersFgmres::GetRequiredMemory(maximumNumberOfNonRestartedIterations, size);
	};
	
      public:
	
	
	virtual std::string GetErrorMessage()const
	{
	  return "Unknown";
	};
	
	//!
	//! @brief Constructor.
	//!
	//! @param size The size of the linear system to solve.
	//! @param matrixVectorProductOperator The matrix vector product operator..
	//! @param numMaxIter The maximum number of iterations.
	//! @param relativeTolerance The relative tolerance.
	//! @param maximumNumberOfNonRestartedIterations The maximum number of non-restarted iterations.
	//! @param hasPreconditioner The matrix vector product operator.
	Fgmres(const wls_int_t  size,
	       ILinearOperator* matrixVectorProductOperator,
	       const wls_int_t  numMaxIter,
	       const double relativeTolerance = 1.0e-6,
	       const wls_int_t  maximumNumberOfNonRestartedIterations = s_maximumNumberOfNonRestartedIterations,
	       const bool hasPreconditioner = true)
	{
       
	  this->m_matrixVectorProductOperator = matrixVectorProductOperator;
	  this->m_preconditioner = NULL;
	  this->m_n = size;
	  this->m_rwork = new double[GetRequiredMemory(size,maximumNumberOfNonRestartedIterations)];

       
	  this->m_status = Status::Completed;
       
	  //
	  // Initialize the solver.
	  //
	  Init(size,
	       &this->m_status,
	       this->m_parameters,
	       this->m_rwork);
       
	  switch (this->m_status)
	    {
	    case Status::Completed:
	      {
		break;
	      }
	    default:
	      {
		std::cerr << "MKLFlexibleGMRES failed:" << Status::GetErrorMessage(m_status) << std::endl;
		exit(1);
	      }
	    }
       
	  //
	  // Personalize the set-up
	  //
       
	  //
	  // Change the default tolerance.
	  //
	  m_parameters.SetRelativeTolerance(relativeTolerance);
       
	  //
	  // Set the maximum number of non-restarted iterations.
	  //
	  m_parameters.SetMaximumNumberOfNonRestartedIterations(maximumNumberOfNonRestartedIterations);
       
	  //
	  // Set the maximum number of iterations. 
	  //
	  m_parameters.SetMaximumNumberOfIterations(numMaxIter);
       
	  //
	  // Let FGMRES do the test for the zero norm of the orthogonal vector.
	  //
	  m_parameters.SetAutomaticTestForZeroNormOfOrthogonalVector(true);
       
	  //
	  // Threshold on the value of the zero norm of the orthogonal vector.
	  //
	  m_parameters.SetZeroNormThresholdOfOrthogonalVector(1.0e-10);
       
	  //
	  // Let FGMRES do the test residual.
	  //
	  m_parameters.SetAutomaticResidualStoppingTest(true);
       
	  //
	  // Let us taking care of the stopping test.
	  //
	  m_parameters.SetAutomaticUserDefinedStoppingTest(false);
       
	  //
	  // Has a preconditioner.
	  //
	  m_parameters.SetPreconditionedMethod(hasPreconditioner);
       
	  Check(this->m_parameters,
		this->m_n,
		this->m_rwork,
		&this->m_status);
       
	  switch (m_status)
	    {
	    case Status::Completed:
	      {
		break;
	      }
	    default:
	      {
		std::cerr << "MKL FGMRES has failed when checking parameters "<<std::endl;
		exit(1);
	      }
	    }
       
	  //
	  // Print information about the set-up.
	  //
	  m_parameters.PrintInformation();
       
	};
	
	//!
	//! @brief Get the parameters.
	//!
	Parameters *GetParameters()
	{
	  return &this->m_parameters;
	};
	
	IInverseOperator* GetPreconditioner()
	{
	  return this->m_preconditioner; 
	};
     
	void SetPreconditioner(IInverseOperator* value)
	{        
	  this->m_preconditioner = value; 
	};

	double m_EuclidianNormOfResidual;
	double GetEuclidianNormOfResidual()const
	{
	  return this->m_EuclidianNormOfResidual;
	};
     
	void SetEuclidianNormOfResidual(const double value)
	{
	  this->m_EuclidianNormOfResidual = value;
	};

	wls_int_t  m_numIterations;
	wls_int_t  GetNumIterations()const
	{
	  return this->m_numIterations;
	};
	void SetNumIterations(const wls_int_t  value)
	{
	  this->m_numIterations = value;
	};

	std::string GetErrorMessage()
	{
	  if (!this->m_preconditionerHasFailed)
	    {
	      this->m_errorMessage = Status::GetErrorMessage(this->m_status);
	    }
	  return this->m_errorMessage;
	};

	void Compute(bool *outHasFailed)
	{
	  *outHasFailed = false;
	  if (this->m_parameters.IsPreconditionedMethod())
	    {
	      if (NULL != this->m_preconditioner)
		{
		  this->m_preconditioner->Compute(&m_preconditionerHasFailed);
		  if (m_preconditionerHasFailed)
		    {
		      *outHasFailed = this->m_preconditionerHasFailed;
		      // this->m_preconditioner->GetErrorMessage;
		      m_errorMessage = std::string("The computation of the preconditioner has failed:");
		    }
		}
	      else
		{
		  std::cerr <<"The preconditioner is missing." << std::endl;
		  exit(1);
		}
	    }
	};
     
     

	wls_int_t  GetSizeOfTemporaryVector() const
	{
	  if (m_parameters.IsPreconditionedMethod())
	    {
	      return m_preconditioner->GetSizeOfTemporaryVector();
	    }
	  else
	    {
	      return 1;
	    }
	};


	void Apply(const char * transpose,
		   double*__restrict__ y,
		   const double*__restrict__ x,
		   const wls_int_t  tmpSize,
		   double*__restrict__ tmp,
		   bool* outHasFailed)
	{
       
	  *outHasFailed = false;
       
	  wls_int_t  numApplyPreconditioner = 0;
	  wls_int_t  numMatrixVectorProduct = 0;
	  //
	  // Checks consistency and correctness of the user defined data. 
	  //
	  this->m_status = Status::Completed;
       
	  //
	  // Makes the FGMRES iterations.
	  //
       
	StateRun:
       
	  Run(this->m_parameters,
	      this->m_n,
	      y,
	      x,
	      this->m_rwork,
	      &this->m_status);
       
	  switch (m_status)
	    {
	    case Status::Completed:
	      {
		break;
	      }
	   
	    case Status::ApplyMatrixVectorProduct:
	      {
		++numMatrixVectorProduct;
		m_matrixVectorProductOperator->Apply(transpose,
						     &this->m_rwork[this->m_parameters.GetIndexOfOutputTemporaryVector()],
						     &this->m_rwork[this->m_parameters.GetIndexOfSourceTemporaryVector()]);
		goto StateRun;
	      }
	   
	    case Status::ApplyPreconditioner:
	      {
#if 0
		Debug.Assert(null != m_preconditioner);
#endif
		++numApplyPreconditioner;
#if 0
		Console.WriteLine("apply preconditioner");
#endif
	      
		m_preconditioner->Apply(transpose,
					&this->m_rwork[this->m_parameters.GetIndexOfOutputTemporaryVector()],
					&this->m_rwork[this->m_parameters.GetIndexOfSourceTemporaryVector()],
					tmpSize,
					tmp,
					&this->m_preconditionerHasFailed);
	     
		if (!m_preconditionerHasFailed)
		  {
		    goto StateRun;
		  }
	     
		*outHasFailed = true;
		m_errorMessage = std::string("The application of the preconditioner has failed:") + m_preconditioner->GetErrorMessage();
		break;
	      }
	   
	    case Status::CheckTheNormOfTheGeneratedVector:
	      {
#if 0
		Debug.Assert(!m_parameters.AutomaticTestForZeroNormOfOrthogonalVector,
			     "The setup of FGMRES is wrong, it says that this step should be automatic.");
#endif
		if (m_parameters.GetNormOfOrthogonalVector() > m_parameters.GetZeroNormThresholdOfOrthogonalVector())
		  {
		    goto StateRun;
		  }

		break;
	      }
	   
	    case Status::MaximumNumberOfIterationsReached:
	      {
		*outHasFailed = true;
		break;
	      }
	   
	    case Status::ApplyStoppingTests:
	      {
	     
		std::cerr << "Apply stopping tests is not yet implemented" << std::endl;
		exit(1);
	      }
	   
	    default:
	      {
		*outHasFailed = true;
		break;
	      }
	    }

	  std::cout
	    << "Status value:" << m_status
	    << ", outHasFailed value: " << *outHasFailed
	    << ", numMatrixVectorProduct: " << numMatrixVectorProduct
	    << ", numApplyPreconditioner: " << numApplyPreconditioner
	    << std::endl;

	
	
#if 0
	  DebugVerbose(string.Format("Status value:{0}, outHasFailed value: {1}", m_status, outHasFailed));
	  Console.WriteLine(
			    "Status value:{0}, outHasFailed value: {1}, numMatrixVectorProduct {2}, numApplyPreconditioner {3}",
			    m_status,
			    outHasFailed,
			    numMatrixVectorProduct,
			    numApplyPreconditioner);
#endif
	  switch (m_status)
	    {
	    case Status::MaximumNumberOfIterationsReached:
	    case Status::Completed:
	      {
	     
		Status::Value requestGet = m_status;
		wls_int_t  itercount = 0;

		Get(m_parameters,
		    m_n,
		    y,
		    (double*__restrict__)x,
		    m_rwork,
		    &requestGet,
		    itercount);
	      
		std::cout << "numIter " << itercount << std::endl;
		

#if 0	     
		Debug.Assert(Status::Completed == requestGet);
#endif	     
		this->SetNumIterations(itercount);
#if 0
		LinearSystemUtils.ComputeResidual(transpose,
						  m_n,
						  y,
						  x,
						  m_rwork,
						  m_matrixVectorProductOperator);
#endif
		static const wls_int_t  n1 = 1;
		this->SetEuclidianNormOfResidual(BlasMKL::nrm2(&m_n, m_rwork,&n1));
#if 0
		DebugVerbose(string.Format("NumIter  = {0}, EuclidianNormOfResidual = {1}", itercount,
					   EuclidianNormOfResidual));
	     
		Console.WriteLine("NumIter  = {0}, EuclidianNormOfResidual = {1}",
				  itercount,
				  EuclidianNormOfResidual);
#endif
	     
		break;
	      }
	   
	    default:
	      {
		this->SetEuclidianNormOfResidual(1e30);
		*outHasFailed = true;
		break;
	      }
	    }
       
	};

     
     

      };
      
    };
    
  };
  
};
