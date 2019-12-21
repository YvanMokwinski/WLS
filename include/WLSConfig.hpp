#pragma once

#if 0
#define _WIN64 0

#if __MNS_ILP__
#define MKL_ILP64 1
#else
#undef MKL_ILP64
#endif
#endif


#define MKL_ILP64 1

#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"

namespace WLS
{
  using integer_t  = long long int;
  using integer_pt = integer_t*__restrict__;
};
