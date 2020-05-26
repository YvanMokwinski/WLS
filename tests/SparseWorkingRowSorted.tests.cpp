

#include "Sparse/WorkingRow/Sorted.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>

class SparseWorkingRowSortedTests : public UnitTest< SparseWorkingRowSortedTests >
{

public:

  inline void run() 
  {    
    this->Test_constructor();
  };
  
  inline void Test_constructor() const 
  {
    static constexpr const WLS::integer_t size = 32;

    WLS::Sparse::WorkingRow::Sorted w(size,
				      size);

#if 0
    Assert::AreEqual(0, w.size());    
    w.add(14);
    Assert::AreEqual(1, w.size());    
#endif    
  };

};

int main()
{
  TestFixture testFixture;
  testFixture.test< SparseWorkingRowSortedTests >();
  return 0;
}


