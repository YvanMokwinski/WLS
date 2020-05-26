#pragma once

namespace WLS
{
  namespace sparse
  {
    
    template <typename _sparsematrix_t,
	      typename _indirect_array_t>
    struct assembly_t
    {
    private: using sparsematrix_t = _sparsematrix_t;
    private: using indirect_array_t = _indirect_array_t;
      
    public: sparsematrix_t& m_sparsematrix;
    public: const indirect_array_t& m_indirect_array;

  public: inline constexpr assembly_t(sparsematrix_t&sparsematrix_,
				    const indirect_array_t& indirect_array_) noexcept
    : m_sparsematrix(sparsematrix_),
      m_indirect_array(indirect_array_)
    {
    };
    
  public: template <typename matrix_t> inline assembly_t& operator += (const matrix_t& that_) noexcept
    {
      const auto size = m_indirect_array.size();
      for (size_t i=0;i<size;++i)
	{
	  const auto idof = m_indirect_array[i];	    
	  for (size_t j=0;j<size;++j)
	    {
	      //	    std::cout << "ass value " << that_(i,j) << std::endl;
	      if (!m_sparsematrix.Ass(idof,m_indirect_array[j],that_(i,j)))
		{
		  std::cerr << "(" << idof << "," <<  m_indirect_array[j] << ") not found" << std::endl;
		  exit(1);
		}
	      //	    if (!m_sparsematrix.Ass(idof,m_indirect_array[j],that_(i,j)))
	      //	      {
	      //		std::cerr << "(" << idof << "," <<  m_indirect_array[j] << ") not found" << std::endl;
	      //	exit(1);
	      // }
	    }
	}    
      return *this;
    };

  };

};
  
};
