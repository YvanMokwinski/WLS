#pragma once

#include "WLS/Sparse/WorkingRow/CRTP.hpp"
namespace WLS
{

  namespace Sparse
  {

    namespace WorkingRow
    {
    
      class Unsorted : public CRTP<Unsorted>
      {
      public: inline Unsorted(const WLS::integer_t size_,
		       const WLS::integer_t subsize_) noexcept;
	
      public: inline ~Unsorted() noexcept;

      public: inline WLS::integer_t size() const noexcept;
      public: inline void add(const WLS::integer_t index_) noexcept;

      public: template <typename _array_t> inline void add(const WLS::integer_t size_,
							   const _array_t& indices_) noexcept;
  
      public: inline void reset() noexcept;
	
      public: inline WLS::integer_t operator[] (const WLS::integer_t localIndex_) const noexcept;

      private: WLS::integer_t 	m_size;
      private: WLS::integer_pt 	m_blank;
      private: WLS::integer_pt 	m_select;
      private: WLS::integer_t 	m_subsize;
  
      };
  
      inline Unsorted::Unsorted(const WLS::integer_t size_,
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
      
      inline Unsorted::~Unsorted() noexcept
      {
	delete [] this->m_blank;
	delete [] this->m_select;
	this->m_subsize = 0;
	this->m_size = 0;
      };
      
      inline WLS::integer_t Unsorted::size() const noexcept
      {
	return this->m_subsize;
      };
      
      inline void Unsorted::add(const WLS::integer_t index_) noexcept
      {
	if (0 == this->m_blank[index_])
	  {
	    this->m_select[this->m_subsize] = index_;
	    this->m_blank[index_] = ++this->m_subsize;
	  }
      };
      
      template <typename _array_t>
      inline void Unsorted::add(const WLS::integer_t size_,
				const _array_t& indices_) noexcept
      {
	for (WLS::integer_t i=0;i<size_;++i)
	  {
	    this->add(indices_[i]);
	  }
      };
      
      inline void Unsorted::reset() noexcept
      {
	for (WLS::integer_t i = 0;i<this->m_subsize;++i)
	  {
	    this->m_blank[this->m_select[i]] = 0;	
	  }
	this->m_subsize = 0;
      };
      
      inline WLS::integer_t Unsorted::operator[] (const WLS::integer_t localIndex_) const noexcept
      {
	return this->m_select[localIndex_];
      };
      
    };
  
  };

};
