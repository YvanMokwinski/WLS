#pragma once


namespace WLS
{
  namespace dense
  {
    struct vector_h;
    struct matrix_h;

    struct temp_gemm
    {
      const double a;
      const char transA;
      const char transB;
      const matrix_h&A;
      const matrix_h&B;
    };
  
    struct temp_gemv
    {
      const double a;
      const char trans;
      const matrix_h&A;
      const vector_h&b;
    };

    struct temp_scal;

    struct temp_transpose
    {
      const char trans;
      const matrix_h&A;
      inline temp_scal operator * (const double&v_) const;
      inline temp_gemv operator * (const vector_h&v_) const;
      inline temp_gemm operator * (const matrix_h&v_) const;
    };

    struct temp_scal
    {
      const double a;
      const char trans;
      const matrix_h&A;
      inline temp_gemv operator * (const vector_h&v_) const;
      inline temp_gemm operator * (const matrix_h&v_) const;
      inline temp_gemm operator * (const temp_transpose&v_) const;
    };

  
    inline temp_gemm temp_scal::operator * (const temp_transpose&v_) const
    {
      return {a,trans,v_.trans,A,v_.A};
    }

    inline temp_gemv temp_scal::operator * (const vector_h&v_) const
    {
      return {a,trans,A, v_};
    };
    inline temp_gemm temp_scal::operator * (const matrix_h&v_) const
    {
      return {a, trans,'N',A, v_};
    };

    inline temp_scal temp_transpose::operator * (const double&v_) const
    {
      return {v_,trans,A};
    };
    inline temp_gemv temp_transpose::operator * (const vector_h&v_) const
    {
      return {1.0,trans,A, v_};
    };
  
    inline temp_gemm temp_transpose::operator * (const matrix_h&v_) const
    {
      return {1.0, trans,'N', A, v_};
    };

    inline temp_scal operator * (const double&v_,const temp_transpose&m) 
    {
      return {v_,m.trans,m.A};
    };
  };
};
