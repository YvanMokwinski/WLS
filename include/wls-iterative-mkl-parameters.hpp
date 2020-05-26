#pragma once

#include "ILinearOperator.hpp"

#include "WLS_MKL.hpp"

namespace WLS
{
  namespace Iterative
  {
    namespace MKL
    {
      
      /// <summary>
      /// Base class for Parameters of  MKL sparse iterative solvers.
      /// </summary>
      class Parameters
      {
	
      protected:
	
	/// <summary>
	/// The array of integer parameters.
	/// </summary>
	wls_int_t m_ipar[128];
	
	/// <summary>
	/// The array of double parameters.
	/// </summary>
	double m_dpar[128];
        
      public:
	
	Parameters()
	{
	  for (int i=0;i<128;++i)
	    {
	      this->m_dpar[i] = ((double)0.0);
	    }
	  for (int i=0;i<128;++i)
	    {
	      this->m_ipar[i] = ((wls_int_t)0);
	    }
	};

      public:
	
	virtual ~Parameters()
	{
	  for (int i=0;i<128;++i)
	    {
	      this->m_dpar[i] = ((double)0.0);
	    }
	  for (int i=0;i<128;++i)
	    {
	      this->m_ipar[i] = ((wls_int_t)0);
	    }
	};
	
      public:
	
	/// <summary>
	/// The array of the integer parameters.
	/// </summary>
	wls_int_p GetParamIntegers()
	{
	  return this->m_ipar;
	};
	
	/// <summary>
	/// The array of the real parameters.
	/// </summary>
	double* GetParamReals()
	{
	  return this->m_dpar;
	};
       
      };
      
    };

  };
  
};
