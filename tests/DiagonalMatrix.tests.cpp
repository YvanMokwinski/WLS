
#include "DiagonalMatrix.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>


class DiagonalMatrixTests : public UnitTest< DiagonalMatrixTests >
{

public:

  inline void run() 
  {    
    this->Test_size<double>();
  };

  template <typename _real_t> inline void Test_size() const 
  {
    WLS::DiagonalMatrix diagonalMatrix(32);
  };
};


int main()
{
  TestFixture testFixture;
  testFixture.test< DiagonalMatrixTests >();
  return 0;
}


