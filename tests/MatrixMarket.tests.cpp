
#include "Direct/MKL/Pardiso.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>
#include <iostream>

#include "WCOMMON/include/cmdline.hpp"
#include "WCOMMON/include/Debug.hpp"
#include "WLA/include/BlasMKL.hpp"

template<typename real_t>
void matrix_market_to_csc(const char * 		filename_,
			  const char * 		idxsys_,
			  WLS::integer_t * 	nrows_,
			  WLS::integer_t * 	ncols_,
			  WLS::integer_t * 	nc_,
			  WLS::integer_t ** 	pbegin_,
			  WLS::integer_t ** 	pidx_,
			  real_t ** 		pvalues_,
			  WLS::integer_t * 	info_)
{  
#ifndef NDEBUG
  Debug::IsNotNull(__DEBUG_TRACE__,filename_);
  Debug::IsNotNull(__DEBUG_TRACE__,idxsys_);
  Debug::IsNotNull(__DEBUG_TRACE__,nrows_);
  Debug::IsNotNull(__DEBUG_TRACE__,ncols_);
  Debug::IsNotNull(__DEBUG_TRACE__,nc_);
  Debug::IsNotNull(__DEBUG_TRACE__,pbegin_);
  Debug::IsNotNull(__DEBUG_TRACE__,pvalues_);
  Debug::IsNotNull(__DEBUG_TRACE__,info_);
#endif
  
  info_[0] = 0;
  {
    FILE * in = fopen(filename_,"r");
    if (!in)
      {
	std::cout << "cannot open file '" << filename_ << "'" << std::endl;
	info_[0] = 1;
	return;
      }
    char*line = nullptr;
    size_t len;
    ssize_t read;
    WLS::integer_t numRows = 0;
    WLS::integer_t numCols = 0;
    WLS::integer_t numCoefficients = 0;
    bool init = false;
    WLS::integer_t numCoefficients_read = 0;

    WLS::integer_t* begin = nullptr;
    WLS::integer_t*idx = nullptr;
    real_t*values = nullptr;

    while ((read = getline(&line, &len, in)) != -1)
      {
	//      printf("Retrieved line of length %zu :\n", read);
	//	printf("%s", line);
	if (line[0] == '%')
	  {
	    //
	    // Comment line.
	    //
	  }
	else
	  {
	    if (!init)
	      {
		init = true;
		sscanf(line,"%Ld %Ld %Ld",&numRows,&numCols,&numCoefficients);
		pbegin_[0] 	= (WLS::integer_t*)calloc((numRows+1),sizeof(WLS::integer_t));
		pidx_[0] 	= (WLS::integer_t*)malloc(sizeof(WLS::integer_t)*(numCoefficients));
		pvalues_[0] = (real_t*)malloc(sizeof(real_t)*(numCoefficients));
		begin = pbegin_[0];
		nrows_[0] = numRows;
		ncols_[0] = numCols;
		nc_[0]    = numCoefficients;
	      }
	    else
	      {
#if 0
		if (!begin)
		  {
		    begin = nullptr;
		  }
#endif
		WLS::integer_t rowIndex = 0;
		WLS::integer_t colIndex = 0;
		double value;
		sscanf(line,"%Ld %Ld %le",&rowIndex,&colIndex,&value);
		if (rowIndex<=0 || colIndex<=0)
		  {
		    info_[0] = 2;
		    fclose(in);
		    return;
		  }
		begin[colIndex-1+1] += 1;
		++numCoefficients_read;
	      }
	  }
    }
    if (line)
      {
	free(line);
	line = nullptr;
      }
    fclose(in);
    in = nullptr;
    
    if (numCoefficients_read != numCoefficients)
      {
	info_[0] = 3;
	return;
      }
    
    begin = pbegin_[0];
    idx = pidx_[0];
    values = pvalues_[0];
		
    //
    //
    //
    for (WLS::integer_t i = 2;i<=numCols;++i)
      {
	begin[i]+=begin[i-1];
      }
    numCoefficients_read = 0;
    
    in = fopen(filename_,"r");
    if (!in)
      {
	std::cout << "cannot open file '" << filename_ << "'" << std::endl;
	info_[0] = 1;
	return;
      }
    
    init = false;
    while ((read = getline(&line, &len, in)) != -1)
      {
	if (line[0] == '%')
	  {
	  }
	else
	  {
	    if (init)
	      {
		WLS::integer_t rowIndex = 0;
		WLS::integer_t colIndex = 0;
		double value;
		sscanf(line,"%Ld %Ld %le\n",&rowIndex,&colIndex,&value);
		if (rowIndex<=0 || colIndex<=0)
		  {
		    info_[0] = 2;
		    fclose(in);
		    return;
		  }
		--colIndex;
		idx[begin[colIndex]] = (idxsys_[0] == 'F') ? rowIndex:rowIndex-1;
		values[begin[colIndex]] = value;
		++begin[colIndex];
		++numCoefficients_read;
	      }
	    else
	      {
		init = true;
	      }
	  }	
    }
    if (line)
      {
	free(line);
	line = nullptr;
      }
    fclose(in);

    if (numCoefficients_read != numCoefficients)
      {
	std::cerr << "error " << numCoefficients_read << " " << numCoefficients << std::endl;
	info_[0] = 3;
	return;
      }

    for (WLS::integer_t i = numCols-1;i>0;--i)
      {
	begin[i]=begin[i-1];
      }
    begin[0]=0;

    //
    // Now we sort the coefficients.
    //
    WLS::integer_t mx = begin[1]-begin[0];
    for (WLS::integer_t i = 1;i<numCols;++i)
      {
	WLS::integer_t n = begin[i+1]-begin[i];
	mx = (mx < n) ? n : mx;
      }

    void * mem = malloc(2*sizeof(double)*mx);
    std::cout << "mx " << mx << " / " << numRows<< std::endl;
    for (WLS::integer_t i = 0;i<numCols;++i)
      {
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    double * m  = (double*)mem;
	    double* x = (double*)(m + 2* (j-begin[i]) + 1);
	    WLS::integer_t* k = (WLS::integer_t*)(x - 1);
	    *k = idx[j];
	    *x = values[j];
	  }
	
	//
	// Can be multi-threaded
	//
	qsort(mem,begin[i+1]-begin[i],2*sizeof(double),
	      [](const void * a,const void * b)
	      {
		WLS::integer_t*a_ = (WLS::integer_t*)a;
		WLS::integer_t*b_ = (WLS::integer_t*)b;
		if (a_[0] > b_[0])
		  {
		    return 1;
		  }
		else if (a_[0] < b_[0])
		  {
		    return -1;
		  }
		else
		  {
		    return 0;
		  }
	      });
	
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    double * m  = (double*)mem;
	    double* x = (double*)(m + 2* (j-begin[i]) + 1);
	    WLS::integer_t* k = (WLS::integer_t*)(x - 1);
	    idx[j] = *k;
	    values[j] = *x;
	  }

#if 0	
	//
	// Can be multi-threaded
	//
	qsort(&idx[begin[i]],begin[i+1]-begin[i],sizeof(WLS::integer_t),
	      [](const void * a,const void * b)
	      {
		WLS::integer_t*a_ = (WLS::integer_t*)a;
		WLS::integer_t*b_ = (WLS::integer_t*)b;
		if (a_[0] < b_[0])
		  {
		    return 1;
		  }
		else if (a_[0] > b_[0])
		  {
		    return -1;
		  }
		else
		  {
		    return 0;
		  }
	      });
#endif	
      }
    free(mem);
    
#if 0
    for (WLS::integer_t i = 0;i<numRows;++i)
      {
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    printf(" %Ld %e",idx[j],values[j]);
	  }
	printf("\n");
      }
#endif
    
    if (idxsys_[0] == 'F')
      {
	for (WLS::integer_t i = 0;i<=numCols;++i)
	  {
	    begin[i]+=1;
	  }
#if 0
	for (WLS::integer_t i = 0;i<numCoefficients;++i)
	  {
	    idx[i]+=1;
	  }
#endif
      }

  }
  
}


template<typename real_t>
void matrix_market_to_csr(const char * 		filename_,
			  const char * 		idxsys_,
			  WLS::integer_t * 	nrows_,
			  WLS::integer_t * 	ncols_,
			  WLS::integer_t * 	nc_,
			  WLS::integer_t ** 	pbegin_,
			  WLS::integer_t ** 	pidx_,
			  real_t ** 		pvalues_,
			  WLS::integer_t * 	info_)
{  
#ifndef NDEBUG
  Debug::IsNotNull(__DEBUG_TRACE__,filename_);
  Debug::IsNotNull(__DEBUG_TRACE__,idxsys_);
  Debug::IsNotNull(__DEBUG_TRACE__,nrows_);
  Debug::IsNotNull(__DEBUG_TRACE__,ncols_);
  Debug::IsNotNull(__DEBUG_TRACE__,nc_);
  Debug::IsNotNull(__DEBUG_TRACE__,pbegin_);
  Debug::IsNotNull(__DEBUG_TRACE__,pvalues_);
  Debug::IsNotNull(__DEBUG_TRACE__,info_);
#endif
  
  info_[0] = 0;
  {
    FILE * in = fopen(filename_,"r");
    if (!in)
      {
	std::cout << "cannot open file '" << filename_ << "'" << std::endl;
	info_[0] = 1;
	return;
      }
    char*line = nullptr;
    size_t len;
    ssize_t read;
    WLS::integer_t numRows = 0;
    WLS::integer_t numCols = 0;
    WLS::integer_t numCoefficients = 0;
    bool init = false;
    WLS::integer_t numCoefficients_read = 0;

    WLS::integer_t* begin = nullptr;
    WLS::integer_t*idx = nullptr;
    real_t*values = nullptr;

    while ((read = getline(&line, &len, in)) != -1)
      {
	//      printf("Retrieved line of length %zu :\n", read);
	//	printf("%s", line);
	if (line[0] == '%')
	  {
	    //
	    // Comment line.
	    //
	  }
	else
	  {
	    if (!init)
	      {
		init = true;
		sscanf(line,"%Ld %Ld %Ld",&numRows,&numCols,&numCoefficients);
		pbegin_[0] 	= (WLS::integer_t*)calloc((numRows+1),sizeof(WLS::integer_t));
		pidx_[0] 	= (WLS::integer_t*)malloc(sizeof(WLS::integer_t)*(numCoefficients));
		pvalues_[0] = (real_t*)malloc(sizeof(real_t)*(numCoefficients));
		begin = pbegin_[0];
		nrows_[0] = numRows;
		ncols_[0] = numCols;
		nc_[0]    = numCoefficients;
		
		//		printf("size %Ld %Ld %Ld\n",numRows,numCols,numCoefficients);
	      }
	    else
	      {
#if 0
		if (!begin)
		  {
		    begin = nullptr;
		  }
#endif
		WLS::integer_t rowIndex = 0;
		WLS::integer_t colIndex = 0;
		double value;
		sscanf(line,"%Ld %Ld %le",&rowIndex,&colIndex,&value);
		if (rowIndex<=0 || colIndex<=0)
		  {
		    info_[0] = 2;
		    fclose(in);
		    return;
		  }
		begin[rowIndex-1+1] += 1;
		++numCoefficients_read;
	      }
	  }
    }
    if (line)
      {
	free(line);
	line = nullptr;
      }
    fclose(in);
    in = nullptr;
    
    if (numCoefficients_read != numCoefficients)
      {
	info_[0] = 3;
	return;
      }
    
    begin = pbegin_[0];
    idx = pidx_[0];
    values = pvalues_[0];
		
    //
    //
    //
    for (WLS::integer_t i = 2;i<=numRows;++i)
      {
	begin[i]+=begin[i-1];
      }
    numCoefficients_read = 0;
    
    in = fopen(filename_,"r");
    if (!in)
      {
	std::cout << "cannot open file '" << filename_ << "'" << std::endl;
	info_[0] = 1;
	return;
      }
    
    init = false;
    while ((read = getline(&line, &len, in)) != -1)
      {
	if (line[0] == '%')
	  {
	  }
	else
	  {
	    if (init)
	      {
		WLS::integer_t rowIndex = 0;
		WLS::integer_t colIndex = 0;
		double value;
		sscanf(line,"%Ld %Ld %le\n",&rowIndex,&colIndex,&value);
		if (rowIndex<=0 || colIndex<=0)
		  {
		    info_[0] = 2;
		    fclose(in);
		    return;
		  }
		--rowIndex;
		idx[begin[rowIndex]] = (idxsys_[0] == 'F') ? colIndex:colIndex-1;
		values[begin[rowIndex]] = value;
		++begin[rowIndex];
		++numCoefficients_read;
	      }
	    else
	      {
		init = true;
	      }
	  }	
    }
    if (line)
      {
	free(line);
	line = nullptr;
      }
    fclose(in);

    if (numCoefficients_read != numCoefficients)
      {
	std::cerr << "error " << numCoefficients_read << " " << numCoefficients << std::endl;
	info_[0] = 3;
	return;
      }

    for (WLS::integer_t i = numRows-1;i>0;--i)
      {
	begin[i]=begin[i-1];
      }
    begin[0]=0;

    //
    // Now we sort the coefficients.
    //
    WLS::integer_t mx = begin[1]-begin[0];
    for (WLS::integer_t i = 1;i<numRows;++i)
      {
	WLS::integer_t n = begin[i+1]-begin[i];
	mx = (mx < n) ? n : mx;
      }

    void * mem = malloc(2*sizeof(double)*mx);
    std::cout << "mx " << mx << " / " << numCols<< std::endl;
    for (WLS::integer_t i = 0;i<numRows;++i)
      {
	BlasMKL::integer_t ln = begin[i+1]-begin[i];
	BlasMKL::integer_t n1=1;
	BlasMKL::integer_t n2=2;
	BlasMKL::copy(&ln,values + begin[i],&n1, ((double*)mem) + 1,&n2);
	
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    WLS::integer_t* k = (WLS::integer_t*)(((double*)mem) + 2* (j-begin[i]) );
	    *k = idx[j];
	  }

#if 0
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  { 
	    double * m  = (double*)mem;
	        double* x = (double*)(m + 2* (j-begin[i]) + 1);
	    WLS::integer_t* k = (WLS::integer_t*)(x - 1);
	    *k = idx[j];
	    *x = values[j];
	  }
#endif	
	//
	// Can be multi-threaded
	//
	qsort(mem,begin[i+1]-begin[i],2*sizeof(double),
	      [](const void * a,const void * b)
	      {
		WLS::integer_t*a_ = (WLS::integer_t*)a;
		WLS::integer_t*b_ = (WLS::integer_t*)b;
		if (a_[0] > b_[0])
		  {
		    return 1;
		  }
		else if (a_[0] < b_[0])
		  {
		    return -1;
		  }
		else
		  {
		    return 0;
		  }
	      });

	BlasMKL::copy(&ln, ((double*)mem) + 1,&n2,values + begin[i],&n1);
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    WLS::integer_t* k = (WLS::integer_t*)(((double*)mem) + 2* (j-begin[i]) );
	    idx[j] = *k;
	  }
	
#if 0
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    double * m  = (double*)mem;
	    double* x = (double*)(m + 2* (j-begin[i]) + 1);
	    WLS::integer_t* k = (WLS::integer_t*)(x - 1);
	    idx[j] = *k;
	    values[j] = *x;
	  }
#endif	
    
#if 0	
	//
	// Can be multi-threaded
	//
	qsort(&idx[begin[i]],begin[i+1]-begin[i],sizeof(WLS::integer_t),
	      [](const void * a,const void * b)
	      {
		WLS::integer_t*a_ = (WLS::integer_t*)a;
		WLS::integer_t*b_ = (WLS::integer_t*)b;
		if (a_[0] < b_[0])
		  {
		    return 1;
		  }
		else if (a_[0] > b_[0])
		  {
		    return -1;
		  }
		else
		  {
		    return 0;
		  }
	      });
#endif	
      }
    free(mem);
    if (idxsys_[0] == 'F')
      {
	for (WLS::integer_t i = 0;i<=numRows;++i)
	  {
	    begin[i]+=1;
	  }
      }

  }
  
}

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
    this->TestSolve<double>();
    this->TestSolve<double>("matrices/e05r0100.mtx");
    
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


  
  template <typename _real_t> inline void TestSolve(const char * ifilename) const 
  {
    WLS::integer_t 	size;
    WLS::integer_t  	ncols;
    WLS::integer_t  	numCoefficients;
    WLS::integer_t * 	begin;
    WLS::integer_t * 	indices;
    double * 		values;
    WLS::integer_t  	info;    
    matrix_market_to_csr<_real_t>(ifilename,
				  "F",
				  &size,
				  &ncols,
				  &numCoefficients,
				  &begin,
				  &indices,
				  &values,
				  &info);
    
    std::cout << "info "  << info << std::endl;
    std::cout << "nrows " << size << std::endl;
    std::cout << "ncols " << ncols << std::endl;
    std::cout << "nc "    << numCoefficients    << std::endl;
    std::cout << "max"    << begin[size] << std::endl;
    
    //
    // Call pardiso.
    //
    WLS::Direct::MKL::Pardiso pardiso(size,
				      begin,
				      indices,
				      values);
    

    //
    // Create expected solution.
    //
    double * expectedSolution = new double[size];
    for (WLS::integer_t i =0;i<size;++i)
      {
	expectedSolution[i] = 1.0;
      }

    //
    // Create RHS.
    //
    double * rhs = new double[size];
    for (WLS::integer_t i = 0;i<size;++i)
      {
	double s =0.0;
	for (WLS::integer_t j = begin[i];j<begin[i+1] ;++j)
	  {
	    s+=values[j-1]*expectedSolution[indices[j-1]-1];
	  }
	rhs[i] = s;
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

#if 0  
    const char * ifilename = cmd.get_arg(1);
  if (verbose)
    {
      std::cout << "ifilename: " << ifilename << std::endl;
    }
  {
    WLS::integer_t 	nrows;
    WLS::integer_t  	ncols;
    WLS::integer_t  	nc;
    WLS::integer_t * 	begin;
    WLS::integer_t * 	idx;
    double * 		values;
    WLS::integer_t  	info;    
    matrix_market_to_csr(ifilename,
			 "C",
			 &nrows,
			 &ncols,
			 &nc,
			 &begin,
			 &idx,
			 &values,
			 &info);
    std::cout << "info "  << info << std::endl;
    std::cout << "nrows " << nrows << std::endl;
    std::cout << "ncols " << ncols << std::endl;
    std::cout << "nc "    << nc    << std::endl;
    std::cout << "max"    << begin[nrows] << std::endl;
  }
#endif  
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


