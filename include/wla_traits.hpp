#pragma once

#include <stdio.h>
#include <iostream>
#include "wls-dense-blas.hpp"
#include <assert.h>

namespace WLS
{
#ifdef WLS_ILP64
  using wla_int_t 		= long long int;
#else
  using wla_int_t 		= int;
#endif
  
  using wla_int_ptr_t 		= wla_int_t *;
  using wla_int_const_ptr_t 	= const wla_int_t *;
};
