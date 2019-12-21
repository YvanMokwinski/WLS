#pragma once
namespace WLS
{
  namespace Sparse
  {

    class SparsityPattern
    {
      
    public: inline bool 		HasFortranIndexing() 	const noexcept;      
    public: inline WLS::integer_t 	GetN() 			const noexcept;
    public: inline WLS::integer_t 	GetM() 			const noexcept;
    public: inline WLS::integer_t 	GetNC() 		const noexcept;
    public: inline WLS::integer_pt 	GetB() 			const noexcept;  
    public: inline WLS::integer_pt 	GetI() 			const noexcept;
      
    public:inline SparsityPattern(const bool 		fortranFormat_,
				  const WLS::integer_t 	n_,
				  const WLS::integer_t 	m_,
				  const WLS::integer_t 	nc_) noexcept;
  
    public:inline SparsityPattern(const bool 		fortranFormat_,
				  const WLS::integer_t 	n_,
				  const WLS::integer_t 	m_,
				  const WLS::integer_t 	nc_,
				  WLS::integer_t*	begin_,
				  WLS::integer_t*	indices_,
				  const bool		owner_) noexcept;

    public: inline WLS::integer_t Find(const WLS::integer_t i_,
				       const WLS::integer_t j_) const noexcept;
#if 0
    public:template <unsigned int _N> static inline void Build(SparsityPattern& sparsityPattern,
							       const std::valarray<std::array<WLS::integer_t,_N> >&m_)
      {
      };
#endif      
      //! @brief number of rows
    private:WLS::integer_t 		m_numRows;
      //! @brief number of columns
    private:WLS::integer_t 		m_numCols;
      //! @brief number of coefficient
    private:WLS::integer_t 		m_numCoefficients;
      //! @brief row bandwith
    private:WLS::integer_t 		m_band_row;
      //! @brief col bandwith
    private:WLS::integer_t 		m_band_col;
      //! @brief format
    private: bool 			m_fortranFormat;
      //! @brief begin index
    private: WLS::integer_t*		m_begin;
      //! @brief indices
    private: WLS::integer_t*		m_indices;
      //! @brief indicate if the owner
    private: bool 			m_owner;

    };


    
    inline bool SparsityPattern::HasFortranIndexing() const noexcept
    {
      return this->m_fortranFormat;
    };

    inline WLS::integer_t SparsityPattern::GetN() const noexcept
    {
      return this->m_numRows;
    };

    inline WLS::integer_t SparsityPattern::GetM() const noexcept
    {
      return this->m_numCols;
    };

    inline WLS::integer_t SparsityPattern::GetNC() const noexcept
    {
      return this->m_numCoefficients;
    };

    inline WLS::integer_pt SparsityPattern::GetB() const noexcept
    {
      return this->m_begin;
    };
  
    inline WLS::integer_pt SparsityPattern::GetI() const noexcept
    {
      return this->m_indices;
    };

    inline WLS::integer_t SparsityPattern::Find(const WLS::integer_t i_,
						const WLS::integer_t j_) const noexcept
    {
      const WLS::integer_t bound = this->m_begin[i_+1];
      if (this->m_fortranFormat)
	{
	  for (WLS::integer_t i=this->m_begin[i_];i<bound;++i)
	    {
	      if (this->m_indices[i-1]==j_+1)
		{
		  return i-1;
		}
	    }
	  return -1;
	}
      else
	{
	  for (WLS::integer_t i=this->m_begin[i_];i<bound;++i)
	    {
	      if (this->m_indices[i]==j_)
		{
		  return i;
		}
	    }
	  return -1;
	}
    };

    inline SparsityPattern::SparsityPattern(const bool 		fortranFormat_,
					    const WLS::integer_t 	n_,
					    const WLS::integer_t 	m_,
					    const WLS::integer_t 	nc_) noexcept
      : m_numRows(n_),
	m_numCols(m_),
	m_numCoefficients(nc_),
	m_fortranFormat(fortranFormat_),
	m_begin(new WLS::integer_t[n_+1]),
	m_indices(new WLS::integer_t[nc_]),
	m_owner(true)
      {
      };
  
    inline SparsityPattern::SparsityPattern(const bool 	fortranFormat_,
					    const WLS::integer_t 	n_,
					    const WLS::integer_t 	m_,
					    const WLS::integer_t 	nc_,
					    WLS::integer_t*	begin_,
					    WLS::integer_t*	indices_,
					    const bool		owner_) noexcept
      : m_numRows(n_),
	m_numCols(m_),
	m_numCoefficients(nc_),
	m_fortranFormat(fortranFormat_),
	m_begin(begin_),
	m_indices(indices_),
	m_owner(owner_)
      {
      };

  };
};
  
