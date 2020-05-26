
#include "Direct/MKL/Pardiso.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>
#include <iostream>
#include "WCOMMON/include/cmdline.hpp"

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
				  double(1.0e-13));
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

    for (int i=0;i<=size_;++i) begin[i]+=1;
    for (int i=0;i<numCoefficients;++i) indices[i]+=1;
    
  };  
};

class PardisoTests : public UnitTest< PardisoTests >
{

public:

  inline void run() 
  {    
    //    this->TestConstructor<double>();
    this->TestSolve<double>();
  };

  template <typename _real_t> inline void TestConstructor() const 
  {
    static constexpr const WLS::integer_t size = 3;
    
    WLS::integer_t 	numCoefficients;
    WLS::integer_pt 	begin;
    WLS::integer_pt 	indices;
    double* 		values;
    
    SparseMatrixGenerator::GenerateTridiagonal(size,
					       &numCoefficients,
					       &begin,
					       &indices,
					       &values);

    for (int i = 0;i < numCoefficients;++i)
      {
	std::cout << values[i] << std::endl;
      }
    
    WLS::Direct::MKL::Pardiso pardiso(size,
				      begin,
				      indices,
				      values);
  };
  
  template <typename _real_t> inline void TestSolve() const 
  {
    static constexpr const WLS::integer_t size = 400;

    WLS::integer_t numCoefficients;
    WLS::integer_pt begin;
    WLS::integer_pt indices;
    double* values;

    //
    // Create matrix.
    //
    SparseMatrixGenerator::GenerateTridiagonal(size,
					       &numCoefficients,
					       &begin,
					       &indices,
					       &values);

    //
    // Call pardiso.
    //
    WLS::Direct::MKL::Pardiso pardiso(size,
				      begin,
				      indices,
				      values);

    

    //
    // Create RHS.
    //
    double * rhs = new double[size];
    rhs[0] = -1.0;
    for (WLS::integer_t i =1;i<size-1;++i)
      {
	rhs[i] = 0.0;
      }
    rhs[size-1] = -1.0;

    //
    // Create expected solution.
    //
    double * expectedSolution = new double[size];
    for (WLS::integer_t i =0;i<size;++i)
      {
	expectedSolution[i] = 1.0;
      }
    
    TestSolveLinearSystem(pardiso,
			  rhs,
			  expectedSolution,
			  size);

  };
  
};


int main(int 		argc,
	 char ** 	argv)
{
  
  WCOMMON::cmdline cmd(argc,argv);

  WLS::integer_t nt;
  if (!cmd.get_integer("--nt",&nt))
    {
      std::cerr << "WARNING, missing '--nt <numthreads>', default is 1." << std::endl;
      nt = 1;
    }

  bool verbose = cmd.get_logical("-v");
  if (verbose)
    {

    }
  
#ifdef WLS_WITH_MKL  
  mkl_set_num_threads(nt);
#endif
  
  TestFixture testFixture;
  testFixture.test< PardisoTests >();
  
#ifdef WLS_WITH_MKL  
  mkl_free_buffers();
  mkl_thread_free_buffers();
#endif  
  return 0;
}


