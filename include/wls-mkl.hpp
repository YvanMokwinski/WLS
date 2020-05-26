#pragma once


#ifdef WLS_WITH_MKL

#define _WIN64 0
#ifdef WLS_ILP64
#define MKL_ILP64 1
#endif

#include "mkl.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "mkl_rci.h"

#endif
