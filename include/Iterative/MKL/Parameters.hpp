#pragma once

#include "ILinearOperator.hpp"

#define MKL_ILP64 1
#include "mkl.h"
#include "mkl_rci.h"


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
	WLS::integer_t m_ipar[128];
	
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
	      this->m_ipar[i] = ((WLS::integer_t)0);
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
	      this->m_ipar[i] = ((WLS::integer_t)0);
	    }
	};
	
      public:
	
	/// <summary>
	/// The array of the integer parameters.
	/// </summary>
	WLS::integer_pt GetParamIntegers()
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
