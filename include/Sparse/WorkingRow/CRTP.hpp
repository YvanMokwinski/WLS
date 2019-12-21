#pragma once

#include "WLSConfig.hpp"

namespace WLS
{
  namespace Sparse
  {
    namespace WorkingRow
    {
      template <typename _impl> class CRTP
      {
      private:   inline const _impl & asImp() const { return static_cast<const _impl&>(*this); }    
      private:   inline  _impl & asImp()  { return static_cast<_impl&>(*this); }  
      public:    inline CRTP(const CRTP<_impl>&) = delete;   
      protected: inline CRTP(){};

      public: inline WLS::integer_t size() const noexcept
	{
	  return asImp().size();
	};

      public: inline void add(const WLS::integer_t index_) noexcept
	{
	  asImp().add(index_);
	};

      public: template <typename _array_t> inline void add(const WLS::integer_t size_,
							   const _array_t& indices_) noexcept
	{
	  asImp().add(size_,
		      indices_);
	};

      public: inline void reset() noexcept
	{
	  asImp().reset();
	};
  
      public: WLS::integer_t operator[] (const WLS::integer_t localIndex_) const noexcept
	{
	  return asImp()[localIndex_];
	};
  
      };
    };
  };
};
