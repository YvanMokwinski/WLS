
#include "WLS/Sparse/WorkingRow/Unsorted.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>

class SparseWorkingRowUnsortedTests : public UnitTest< SparseWorkingRowUnsortedTests >
{

public:

  inline void run() 
  {    
    this->Test_size<double>();
  };
  
  template <typename _real_t> inline void Test_size() const 
  {
    WLS::Sparse::WorkingRow::Unsorted w(32,32);
  };
};

int main()
{
  TestFixture testFixture;
  testFixture.test< SparseWorkingRowUnsortedTests >();
  return 0;
}


