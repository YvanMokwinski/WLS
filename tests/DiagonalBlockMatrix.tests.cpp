
#include "DiagonalBlockMatrix.hpp"
#include "Assert.hpp"
#include "TestFixture.hpp"
#include <array>




class DiagonalBlockMatrixTests : public UnitTest< DiagonalBlockMatrixTests >
{

public:

  inline void run() 
  {    
    this->Test_size<double>();
  };

  
  template <typename _real_t> inline void Test_size() const 
  {
    WLA::matrix block(2,2);
    
    static_cast<WLA::matrix_h&>(block) = 0.0; 
    block.x[0] = 1.0;
    block.x[3] = 2.0;
    WLS::integer_t size_block = 2;
    WLS::integer_t nblocks    = 32;
    WLS::DiagonalBlockMatrixFactory diagonalBlockMatrixFactory(nblocks,size_block);
    for (WLS::integer_t block_index = 0;block_index < nblocks;++block_index)
      {
	diagonalBlockMatrixFactory.addelm(block_index,1.0, block);	
      }       


    
    
    WLS::DiagonalBlockMatrix * inverseDiagonalBlockMatrix = diagonalBlockMatrixFactory.create_inverse();
    std::cout << *inverseDiagonalBlockMatrix << std::endl;

    delete inverseDiagonalBlockMatrix;
  };

  
};


int main()
{
  TestFixture testFixture;
  testFixture.test< DiagonalBlockMatrixTests >();
  return 0;
}


