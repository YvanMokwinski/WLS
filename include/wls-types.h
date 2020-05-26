#ifndef WLS_TYPES_H
#define WLS_TYPES_H

#ifndef WLS_ILP64
#error WLS_ILP64 is undefined.
#endif

//!
//! @brief Define integer type.
//!
#if WLS_ILP64

typedef long long int wls_int_t;
#define WLS_INT_FORMAT "%Ld"

#else

typedef int wls_int_t;
#define WLS_INT_FORMAT "%d"

#endif

//!
//! @brief Define integer pointer type.
//!
typedef wls_int_t * __restrict__ 	wls_int_p;

//!
//! @brief Define unmutable integer pointer type.
//!
typedef const wls_int_t * __restrict__ const_wls_int_p;


//!
//! @brief Define string type.
//!
typedef char wls_str_t[256];

#ifdef WLS_ILP64
typedef long long int integer_t;
#else
typedef int integer_t;
#endif
typedef wls_int_t wls_status_t;



#define iformat "%ld"
#ifdef WLS_ILP64
#undef  iformat
#define iformat "%Ld"
#endif

#endif

