#pragma once

namespace WLS
{
  namespace Sparse
  {
    /// <summary>
    /// Routines for the SparseBLAS.
    /// </summary>
    class SpBlas
    {
    public: using integer_t = long long int;
    public: using integer_pt = long long int*;
    
      /// <summary>
      /// General matrix vector product for a sparse matrix with the Compressed Sparse Row format with 1-based indexing.
      /// </summary>
      /// <param name="transa">No transpose ('N' or 'n'), Transpose ('T' or 't')</param>
      /// <param name="m">Number of rows of the matrix A.</param>
      /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
      /// <param name="ia">Array of length <paramref name="m"/> + 1, containing indices of elements in the array a, such that <paramref name="ia"/>[I] is the index in the array a of the first non-zero element from the row I. The value of the last element <paramref name="ia"/>[m] is equal to the number of non-zeros.</param>
      /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A, its length is equal to the length of the array <paramref name="a"/>. </param>
      /// <param name="x">Input vector.</param>
      /// <param name="y">Output vector.</param>
    public: template <typename real_t> static void csrgemv1based(const char* transa_,
								 const integer_t m_,
								 const real_t * a_,
								 const integer_pt ia_,
								 const integer_pt ja_,
								 const real_t*x_,
								 real_t*y_);
    
      /// <summary>
      /// Triangular solvers with simplified interface for a sparse matrix in the CSR format (3-array variation) with one-based indexing.
      /// </summary>
      /// <param name="uplo">Specifies whether the upper or low triangle of the matrix A is used. 
      /// If <paramref name="uplo"/> = 'U' or 'u', then the upper triangle of the matrix A is used. 
      /// If <paramref name="uplo"/> = 'L' or 'l', then the low triangle of the matrix A is used.</param>
      /// <param name="transa">Specifies the system of linear equations.
      /// If <paramref name="transa"/> = 'N' or 'n', then A*y = x.
      /// If <paramref name="transa"/> = 'T' or 't' or 'C' or 'c', then AT*y = x,</param>
      /// <param name="diag">Specifies whether A is unit triangular.
      /// If <paramref name="diag"/> = 'U' or 'u', then A is a unit triangular. 
      /// If <paramref name="diag"/> = 'N' or 'n', then A is not unit triangular.</param>
      /// <param name="m"> Number of rows of the matrix A.</param>
      /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
      /// <param name="ia">Array of length If <paramref name="m"/> + 1, containing indices of elements in the array <paramref name="a"/>, 
      /// such that <paramref name="ia"/>(i) is the index in the array <paramref name="a"/> of the first non-zero element from the row i. 
      /// The value of the last element <paramref name="ia"/>(<paramref name="m"/> + 1) is equal to the number of non-zeros plus one. </param>
      /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array a.</param>
      /// <param name="x">Array, size is <paramref name="m"/>. On entry, the array x must contain the vector x.</param>
      /// <param name="y">Array, size at least <paramref name="m"/>. Contains the vector y.</param>
      /// <remarks>
      /// 
      /// The mkl_dcsrtrsv routine solves a system of linear equations with matrix-vector operations for a sparse matrix stored in the CSR format (3 array variation):
      /// A*<paramref name="y"/> = <paramref name="x"/> or A^T * <paramref name="y"/> = <paramref name="x"/>.
      /// where: <paramref name="x"/> and <paramref name="y"/> are vectors, A is a sparse upper or lower triangular matrix with unit or non-unit main diagonal, A^T is the transpose of A.
      /// 
      /// IMPORTANT: This routine supports only one-based indexing of the input arrays.
      /// </remarks>
    public: template <typename real_t> static void csrtrsv1based(const char* uplo_,
								 const char* transa_,
								 const char* diag_,
								 const integer_t m_,
								 const real_t* a_,
								 const integer_pt ia_,
								 const integer_pt ja_,
								 const real_t* x_,
								 real_t* y_);
    


      /// <summary>
      /// Triangular solvers with simplified interface for a sparse matrix in the CSR format (3-array variation) with zero-based indexing.
      /// </summary>
      /// <param name="uplo">Specifies whether the upper or low triangle of the matrix A is used. 
      /// If <paramref name="uplo"/> = 'U' or 'u', then the upper triangle of the matrix A is used. 
      /// If <paramref name="uplo"/> = 'L' or 'l', then the low triangle of the matrix A is used.</param>
      /// <param name="transa">Specifies the system of linear equations.
      /// If <paramref name="transa"/> = 'N' or 'n', then A*y = x.
      /// If <paramref name="transa"/> = 'T' or 't' or 'C' or 'c', then AT*y = x,</param>
      /// <param name="diag">Specifies whether A is unit triangular.
      /// If <paramref name="diag"/> = 'U' or 'u', then A is a unit triangular. 
      /// If <paramref name="diag"/> = 'N' or 'n', then A is not unit triangular.</param>
      /// <param name="m"> Number of rows of the matrix A.</param>
      /// <param name="a">Array containing non-zero elements of the matrix A. Its length is equal to the number of non-zero elements in the matrix A.</param>
      /// <param name="ia">Array of length If <paramref name="m"/> + 1, containing indices of elements in the array <paramref name="a"/>, 
      /// such that <paramref name="ia"/>(i) is the index in the array <paramref name="a"/> of the first non-zero element from the row i. 
      /// The value of the last element <paramref name="ia"/>(<paramref name="m"/> + 1) is equal to the number of non-zeros plus one. </param>
      /// <param name="ja">Array containing the column indices for each non-zero element of the matrix A. Its length is equal to the length of the array a.</param>
      /// <param name="x">Array, size is <paramref name="m"/>. On entry, the array x must contain the vector x.</param>
      /// <param name="y">Array, size at least <paramref name="m"/>. Contains the vector y.</param>
      /// <remarks>
      /// 
      /// The mkl_dcsrtrsv routine solves a system of linear equations with matrix-vector operations for a sparse matrix stored in the CSR format (3 array variation):
      /// A*<paramref name="y"/> = <paramref name="x"/> or A^T * <paramref name="y"/> = <paramref name="x"/>.
      /// where: <paramref name="x"/> and <paramref name="y"/> are vectors, A is a sparse upper or lower triangular matrix with unit or non-unit main diagonal, A^T is the transpose of A.
      /// 
      /// IMPORTANT: This routine supports only one-based indexing of the input arrays.
      /// </remarks>
    public: template <typename real_t> static void csrtrsv0based(const char* uplo_,
								 const char* transa_,
								 const char* diag_,
								 const integer_t m_,
								 const real_t* a_,
								 const integer_pt ia_,
								 const integer_pt ja_,
								 const real_t* x_,
								 real_t* y_);
    
    };
  
    template <> void SpBlas::csrgemv1based<double>(const char* transa_,
						   const SpBlas::integer_t m_,
						   const double * a_,
						   const SpBlas::integer_pt ia_,
						   const SpBlas::integer_pt ja_,
						   const double*x_,
						   double*y_)
    {
      mkl_dcsrgemv((char*)transa_,
		   (SpBlas::integer_pt)&m_,
		   (double*)a_,
		   (SpBlas::integer_pt)ia_,
		   (SpBlas::integer_pt)ja_,
		   (double*)x_,
		   (double*)y_);
    };

    template <> void SpBlas::csrgemv1based<float>(const char* transa_,
						  const SpBlas::integer_t m_,
						  const float * a_,
						  const SpBlas::integer_pt ia_,
						  const SpBlas::integer_pt ja_,
						  const float*x_,
						  float*y_)
    {
#if 0
      mkl_scsrgemv((char*)transa_,
		   (SpBlas::integer_t)&m_,
		   (float*)a_,
		   (SpBlas::integer_pt)ia_,
		   (SpBlas::integer_pt)ja_,
		   (float*)x_,
		   (float*)y_);
#endif
      std::cerr << "come here" << std::endl;
    };



    template <> void SpBlas::csrtrsv0based<double>(const char* uplo_,
						   const char* transa_,
						   const char* diag_,
						   const SpBlas::integer_t m_,
						   const double* a_,
						   const SpBlas::integer_pt ia_,
						   const SpBlas::integer_pt ja_,
						   const double* x_,
						   double* y_)
    {
      mkl_cspblas_dcsrtrsv((char*)uplo_,
			   (char*)transa_,
			   (char*)diag_,
			   (SpBlas::integer_pt)&m_,
			   (double*)a_,
			   (SpBlas::integer_pt) ia_,
			   (SpBlas::integer_pt) ja_,
			   (double*) x_,
			   (double*) y_);    
    };
  
    template <> void SpBlas::csrtrsv0based<float>(const char* uplo_,
						  const char* transa_,
						  const char* diag_,
						  const SpBlas::integer_t m_,
						  const float* a_,
						  const SpBlas::integer_pt ia_,
						  const SpBlas::integer_pt ja_,
						  const float* x_,
						  float* y_)
    {
      mkl_cspblas_scsrtrsv((char*)uplo_,
			   (char*)transa_,
			   (char*)diag_,
			   (SpBlas::integer_pt)&m_,
			   (float*)a_,
			   (SpBlas::integer_pt) ia_,
			   (SpBlas::integer_pt) ja_,
			   (float*) x_,
			   (float*) y_);    
    };


    template <> void SpBlas::csrtrsv1based<double>(const char* uplo_,
						   const char* transa_,
						   const char* diag_,
						   const SpBlas::integer_t m_,
						   const double* a_,
						   const SpBlas::integer_pt ia_,
						   const SpBlas::integer_pt ja_,
						   const double* x_,
						   double* y_)
    {
      mkl_dcsrtrsv((char*)uplo_,
		   (char*)transa_,
		   (char*)diag_,
		   (SpBlas::integer_pt)&m_,
		   (double*)a_,
		   (SpBlas::integer_pt) ia_,
		   (SpBlas::integer_pt) ja_,
		   (double*) x_,
		   (double*) y_);    
    };

    template <> void SpBlas::csrtrsv1based<float>(const char* uplo_,
						  const char* transa_,
						  const char* diag_,
						  const SpBlas::integer_t m_,
						  const float* a_,
						  const SpBlas::integer_pt ia_,
						  const SpBlas::integer_pt ja_,
						  const float* x_,
						  float* y_)
    {
      mkl_scsrtrsv((char*)uplo_,
		   (char*)transa_,
		   (char*)diag_,
		   (SpBlas::integer_pt)&m_,
		   (float*)a_,
		   (SpBlas::integer_pt) ia_,
		   (SpBlas::integer_pt) ja_,
		   (float*) x_,
		   (float*) y_);    
    };

  };

};
