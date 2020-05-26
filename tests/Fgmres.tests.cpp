#include <iostream>

#include "Assert.hpp"
#include "TestFixture.hpp"

#include "wls-direct-mkl-pardiso.hpp"
#include "wls-iterative-mkl-fgmres.hpp"
#include "wls-sparse-matrix.hpp"
#include "wls-preconditioner-mkl-ilu0.hpp"
#include "wls-preconditioner-mkl-ilut.hpp"

void TestSolveLinearSystem(WLS::IInverseOperator& 	inverseOperator_,
			   const double * 		rhs_,
			   const double * 		expectedSolution_,
			   const WLS::integer_t 	size_)
{
  
  double * solution = new double[size_];
  
  bool hasFailed;
  inverseOperator_.Compute(&hasFailed);
  if (hasFailed)
    {
      std::cerr << "compute failed: " << inverseOperator_.GetErrorMessage()  << std::endl;
      exit(1);
    }
  
  double* temporaryVector = nullptr;
  WLS::integer_t sizeOfTemporaryVector = inverseOperator_.GetSizeOfTemporaryVector();
  if (sizeOfTemporaryVector>0)
    {
      temporaryVector = new double[sizeOfTemporaryVector];
    }

  
  inverseOperator_.Apply("No transpose",
			 solution,
			 rhs_,
			 sizeOfTemporaryVector,
			 temporaryVector,
			 &hasFailed);

  for (int i=0;i<size_;++i)
    {
      Assert::AreNumericallyEqual(expectedSolution_[i],
				  solution[i],
				  double(0.000000001));
    }
  
  if (hasFailed)
    {
      std::cerr << "apply failed: " << inverseOperator_.GetErrorMessage()  << std::endl;
      exit(1);
    }

  if (nullptr != solution)
    {
      delete[] solution;
      solution = nullptr;
    }
  
  if (nullptr != temporaryVector)
    {
      delete[] temporaryVector;
      temporaryVector = nullptr;
    }

};


class SparseMatrixGenerator
{
  
public: static void GenerateTridiagonal(const WLS::integer_t  	size_,
					WLS::integer_t * 	outNumCoefficients_,
					WLS::integer_pt * 	outBegin_,
					WLS::integer_pt * 	outIndices_,
					double ** 		outValues_) noexcept
  {    
    WLS::integer_t 	numCoefficients = 3*(size_-2) + 2*2;
    WLS::integer_pt 	begin = new WLS::integer_t[size_+1];
    const WLS::integer_t n = size_-1;
    
    {
      begin[0] = 0;
      begin[1] = 2;
      for (WLS::integer_t i = 2;i < size_;++i)
	{
	  begin[i] = begin[i-1] + 3;
	}
      begin[size_] = begin[size_-1] + 2;
    }

    WLS::integer_pt indices = new WLS::integer_t[numCoefficients];
    {
      indices[0] = 0;
      indices[1] = 1;
      WLS::integer_t at = 2;    
      for (WLS::integer_t i = 1; i <n; ++i)
	{
	  indices[at++] = i - 1;
	  indices[at++] = i;
	  indices[at++] = i + 1;
	}
      indices[at++] = n - 1;
      indices[at++] = n;      
    }

    double* values = new double[numCoefficients];
    {
      static constexpr const double one = double(1.0);
      static constexpr const double mtwo = double(-2.0);

      values[0] = mtwo;
      values[1] = one;
      WLS::integer_t at = 2;    
      for (WLS::integer_t i = 1; i <n; ++i)
	{
	  values[at++] = one;
	  values[at++] = mtwo;
	  values[at++] = one;
	}
      values[at++] = one;
      values[at++] = mtwo;      
    }
    
    outNumCoefficients_[0] = numCoefficients;
    outBegin_[0] = begin;
    outIndices_[0] = indices;    
    outValues_[0] = values;

    for (WLS::integer_t i=0;i<=size_;++i) begin[i]+=1;
    for (WLS::integer_t i=0;i<numCoefficients;++i) indices[i]+=1;

  };  
};

class FgmresTests : public UnitTest< FgmresTests >
{

public:

  inline void run() 
  {
    this->TestSolve<double>();
  };

  
  template <typename _real_t> inline void TestSolve() const 
  {
    
    static constexpr const WLS::integer_t size = 3000;

    WLS::integer_t numCoefficients;
    WLS::integer_pt begin;
    WLS::integer_pt indices;
    double* values;
    
    SparseMatrixGenerator::GenerateTridiagonal(size,
					       &numCoefficients,
					       &begin,
					       &indices,
					       &values);


    WLS::Sparse::SparsityPattern* sparsityPattern = new WLS::Sparse::SparsityPattern(true,
										     size,size,
										     numCoefficients,
										     begin,
										     indices,
										     true);
    WLS::Sparse::Matrix<double> matrix(sparsityPattern);
    double*x = matrix.GetX();
    for (WLS::integer_t i=0;i<numCoefficients;++i)
      {
	x[i] = values[i];
      }
    
    WLS::Iterative::MKL::Fgmres linearSolver(size,
					     matrix.GetLinearOperator(),
					     10000,
					     1.0e-13,
					     30,
					     true);
    
#if 1
    WLS::Preconditioner::MKL::Ilu0 ilu0(&matrix,
					linearSolver.GetParameters());
    linearSolver.SetPreconditioner(&ilu0);
#endif
#if 0
    WLS::Preconditioner::MKL::Ilut ilut(&matrix,
					0.01,
					3,
					linearSolver.GetParameters());
    
    
    // 
    linearSolver.SetPreconditioner(&ilut);
#endif
    
    double * rhs = new double[size];
    rhs[0] = -1.0;
    for (WLS::integer_t i =1;i<size-1;++i)
      {
	rhs[i] = 0.0;
      }
    rhs[size-1] = -1.0;

    double * expectedSolution = new double[size];
    for (WLS::integer_t i =0;i<size;++i)
      {
	expectedSolution[i] = 1.0;
      }

    TestSolveLinearSystem(linearSolver,
			  rhs,
			  expectedSolution,
			  size);

  };

  
};


int main()
{
#ifdef WLS_WITH_MKL
  mkl_set_num_threads(4);
#endif
  TestFixture testFixture;
  testFixture.test< FgmresTests >();
#ifdef WLS_WITH_MKL  
  mkl_free_buffers();
  mkl_thread_free_buffers();
#endif
  return 0;
}


