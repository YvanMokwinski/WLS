#pragma once
namespace WLS
{
  namespace sparse
  {

    class symbolic_t
    {
      
    public: inline bool 		HasFortranIndexing() 	const noexcept;      
    public: inline wls_int_t 	GetN() 			const noexcept;
    public: inline wls_int_t 	GetM() 			const noexcept;
    public: inline wls_int_t 	GetNC() 		const noexcept;
    public: inline wls_int_p 	GetB() 			const noexcept;  
    public: inline wls_int_p 	GetI() 			const noexcept;
      
    public: inline symbolic_t(const bool 			fortranFormat_,
				   wls_int_t 	n_,
				   wls_int_t 	m_,
				   wls_int_t 	nc_) noexcept;
      
    public: inline symbolic_t(const bool 			fortranFormat_,
				   wls_int_t 	n_,
				   wls_int_t 	m_,
				   wls_int_t 	nc_,
				   wls_int_t*		begin_,
				   wls_int_t*		indices_,
				   const bool			owner_) noexcept;
      
    public: inline wls_int_t Find(wls_int_t i_,
				       wls_int_t j_) const noexcept;
#if 0
    public:template <unsigned int _N> static inline void Build(symbolic_t& sparsityPattern,
							       const std::valarray<std::array<wls_int_t,_N> >&m_)
      {
      };
#endif      
      //! @brief number of rows
    private:wls_int_t 		m_numRows;
      //! @brief number of columns
    private:wls_int_t 		m_numCols;
      //! @brief number of coefficient
    private:wls_int_t 		m_numCoefficients;
      //! @brief row bandwith
    private:wls_int_t 		m_band_row;
      //! @brief col bandwith
    private:wls_int_t 		m_band_col;
      //! @brief format
    private: bool 			m_fortranFormat;
      //! @brief begin index
    private: wls_int_t*		m_begin;
      //! @brief indices
    private: wls_int_t*		m_indices;
      //! @brief indicate if the owner
    private: bool 			m_owner;

    };


    
    inline bool symbolic_t::HasFortranIndexing() const noexcept
    {
      return this->m_fortranFormat;
    };

    inline wls_int_t symbolic_t::GetN() const noexcept
    {
      return this->m_numRows;
    };

    inline wls_int_t symbolic_t::GetM() const noexcept
    {
      return this->m_numCols;
    };

    inline wls_int_t symbolic_t::GetNC() const noexcept
    {
      return this->m_numCoefficients;
    };

    inline wls_int_p symbolic_t::GetB() const noexcept
    {
      return this->m_begin;
    };
  
    inline wls_int_p symbolic_t::GetI() const noexcept
    {
      return this->m_indices;
    };

    inline wls_int_t symbolic_t::Find(wls_int_t i_,
						wls_int_t j_) const noexcept
    {
      wls_int_t bound = this->m_begin[i_+1];
      if (this->m_fortranFormat)
	{
	  for (wls_int_t i=this->m_begin[i_];i<bound;++i)
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
	  for (wls_int_t i=this->m_begin[i_];i<bound;++i)
	    {
	      if (this->m_indices[i]==j_)
		{
		  return i;
		}
	    }
	  return -1;
	}
    };

    inline symbolic_t::symbolic_t(const bool 		fortranFormat_,
					    wls_int_t 	n_,
					    wls_int_t 	m_,
					    wls_int_t 	nc_) noexcept
      : m_numRows(n_),
	m_numCols(m_),
	m_numCoefficients(nc_),
	m_fortranFormat(fortranFormat_),
	m_begin(new wls_int_t[n_+1]),
	m_indices(new wls_int_t[nc_]),
	m_owner(true)
      {
      };
  
    inline symbolic_t::symbolic_t(const bool 	fortranFormat_,
					    wls_int_t 	n_,
					    wls_int_t 	m_,
					    wls_int_t 	nc_,
					    wls_int_t*	begin_,
					    wls_int_t*	indices_,
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
  
