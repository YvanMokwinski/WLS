#pragma once

#include "WLS/Sparse/WorkingRow/CRTP.hpp"
namespace WLS
{

  namespace Sparse
  {

    namespace WorkingRow
    {
      class Sorted : public CRTP<Sorted>
      {
      private: WLS::integer_t 	m_size;
      private: WLS::integer_pt 	m_blank;
      private: WLS::integer_pt 	m_select;
      private: WLS::integer_t 	m_subsize;
  
      public: inline Sorted(const WLS::integer_t size_,
		     const WLS::integer_t subsize_) noexcept;

      public: inline ~Sorted() noexcept;

      public: inline WLS::integer_t size() const noexcept;
	
      public: inline void add(const WLS::integer_t index_) noexcept;
	
      public: template <typename _array_t>
      inline void add(const WLS::integer_t size_,
		      const _array_t& indices_) noexcept;
	
      public: inline void reset() noexcept;
  
      public: inline WLS::integer_t operator[] (const WLS::integer_t localIndex_) const noexcept;  
      };



  
      inline Sorted::Sorted(const WLS::integer_t size_,
			    const WLS::integer_t subsize_) noexcept
	: m_size(size_),
	  m_blank(new WLS::integer_t[size_]),
	  m_select(new WLS::integer_t[subsize_]),
	  m_subsize(0)
      {
	for (WLS::integer_t i=0;i<size_;++i)
	  {
	    m_blank[i] = 0;
	  }
      };
      
      inline Sorted::~Sorted() noexcept
      {
	delete [] this->m_blank;
	delete [] this->m_select;
	this->m_subsize = 0;
	this->m_size = 0;
      };
      
      inline WLS::integer_t Sorted::size() const noexcept
      {
	return this->m_subsize;
      };
      
      inline void Sorted::add(const WLS::integer_t index_) noexcept
      {
	if (0 == this->m_blank[index_])
	  {
	    if (m_subsize>0)
	      {	    
		bool found = false;
		if (index_ < this->m_select[this->m_subsize-1])
		  {
		    for (WLS::integer_t i=0;i<this->m_subsize;++i)
		      {
			if (index_ < this->m_select[i])
			  {
			    for (WLS::integer_t l=this->m_subsize;l>i;--l)
			      {
				this->m_select[l] = this->m_select[l-1];
			      }
			    found = true;
			    this->m_select[i] = index_;		  
			    break;
			  }
		      }
		  }
		
		if (!found)
		  {
		    this->m_select[this->m_subsize] = index_;
		  }
	      }
	    else
	      {
		this->m_select[this->m_subsize] = index_; 
	      }
	    this->m_blank[index_] = ++this->m_subsize;
	  }
      };
      
      template <typename _array_t>
      inline void Sorted::add(const WLS::integer_t size_,
			      const _array_t& indices_) noexcept
      {
	for (WLS::integer_t i=0;i<size_;++i)
	  {
	    this->add(indices_[i]);
	  }
      };
      
      inline void Sorted::reset() noexcept
      {
	for (WLS::integer_t i = 0;i<this->m_subsize;++i)
	  {
	    this->m_blank[this->m_select[i]] = 0;	
	  }
	this->m_subsize = 0;
      };
      
      inline WLS::integer_t Sorted::operator[] (const WLS::integer_t localIndex_) const noexcept
      {
	return this->m_select[localIndex_];
      };
      
    };
    
  };

};