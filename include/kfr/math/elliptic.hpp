/** @addtogroup other_math
 *  @{
 */

#pragma once

#include "impl/elliptic.hpp"
#include <iostream>
#include <complex>

namespace kfr
{
inline namespace CMT_ARCH_NAME
{

//****************************************************************************80

/* 
 * Some method used for conversion and calculation
 * @TODO To be rewriten maybe.. kfr style ? 
 */
template <typename T>
complex<T> sqrt(complex<T> w) 
{
  std::complex<T> z = std::sqrt(std::complex(w.real(), w.imag()));
  return kfr::complex(z.real(), z.imag());
}
template <typename T>
complex<T> atanh(complex<T> w) 
{
  std::complex<T> z = std::atanh(std::complex(w.real(), w.imag()));
  return kfr::complex(z.real(), z.imag());
}
template <typename T>
complex<T> asin(complex<T> w) 
{
  std::complex<T> z = std::asin(std::complex(w.real(), w.imag()));
  return kfr::complex(z.real(), z.imag());
}

template <typename T>
T cmpl(T kx) { return sqrt((1 - kx) * (1 + kx)); }

template <typename T>
complex<T> cmpl(complex<T> kx) { return sqrt((1 - kx) * (1 + kx)); }

template <typename T>
T cephes_polevl(T x, T *coef, int N)
{
    T	ans;
    int		i;
    T	*p;

    p = coef;
    ans = *p++;
    i = N;

    do
      ans = ans * x  +  *p++;
    while ( --i );

    return ans;
}


//****************************************************************************80

template <typename T>
void sncndn ( T u, T m, T &sn, T &cn, T &dn )

//****************************************************************************80
//
//  Purpose:
//
//    SNCNDN evaluates Jacobi elliptic functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    Original ALGOL version by Roland Bulirsch.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Roland Bulirsch,
//    Numerical calculation of elliptic integrals and elliptic functions,
//    Numerische Mathematik,
//    Volume 7, Number 1, 1965, pages 78-90.
//
//  Parameters:
//
//    Input, T U, M, the arguments.
//
//    Output, T &SN, &CN, &DN, the value of the Jacobi
//    elliptic functions sn(u,m), cn(u,m), and dn(u,m).
//
{
  T a;
  T b;
  T c;
  T ca;
  T d;
  int i;
  int l;
  T *m_array;
  T m_comp;
  T *n_array;
  const T r8_epsilon = 2.220446049250313E-16;
  T u_copy;

  m_comp = 1.0 - m;
  u_copy = u;

  if ( m_comp == 0.0 )
  {
    cn = 1.0 / cosh ( u_copy );
    dn = cn;
    sn = tanh ( u_copy );
    return;
  }

  if ( 1.0 < m )
  {
    d = 1.0 - m_comp;
    m_comp = - m_comp / d;
    d = sqrt ( d );
    u_copy = d * u_copy;
  }

  ca = sqrt ( r8_epsilon );

  a = 1.0;
  dn = 1.0;
  l = 24;

  m_array = new double[25];
  n_array = new double[25];

  for ( i = 0; i < 25; i++ )
  {
    m_array[i] = a;
    m_comp = sqrt ( m_comp );
    n_array[i] = m_comp;
    c = 0.5 * ( a + m_comp );
    if ( fabs ( a - m_comp ) <= ca * a )
    {
      l = i;
      break;
    }
    m_comp = a * m_comp;
    a = c;
  }

  u_copy = c * u_copy;
  sn = sin ( u_copy );
  cn = cos ( u_copy );

  if ( sn != 0.0 )
  {
    a = cn / sn;
    c = a * c;

    for ( i = l; 0 <= i; i-- )
    {
      b = m_array[i];
      a = c * a;
      c = dn * c;
      dn = ( n_array[i] + a ) / ( b + a );
      a = c / b;
    }

    a = 1.0 / sqrt ( c * c + 1.0 );

    if ( sn < 0.0 )
    {
      sn = - a;
    }
    else
    {
      sn = a;
    }
    cn = c * sn;
  }

  if ( 1.0 < m )
  {
    a = dn;
    dn = cn;
    cn = a;
    sn = sn / d;
  }

  delete [] m_array;
  delete [] n_array;

  return;
}

//****************************************************************************80

/**  Inverse Jacobian elliptic sn

    See [1], Eq. (56)

    References
    ----------
    .. [1] Orfanidis, "Lecture Notes on Elliptic Filter Design",
           https://www.ece.rutgers.edu/~orfanidi/ece521/notes.pdf

*/

template <typename T>
complex<T> inv_jacobi_sn(complex<T> w, T m) 
{
    /* Maximum number of iterations in Landen transformation recursion
       sequence.  10 is conservative; unit tests pass with 4, Orfanidis
       (see jacobi_cn [1]) suggests 5. */
    const int INV_JACOBI_SN_MAXITER = 10;

    T k = sqrt(m);
    if(k > 1) return NAN;
    else if( k == 1) return atanh(w);

    int niter = 0;
    univector<T> ks(1, k);
    while (ks[ks.size() - 1] != 0)
    {
        T k_  = ks[ks.size() - 1];
        T k_p = cmpl(k_);
        ks.push_back((1 - k_p) / (1 + k_p));
  
        niter += 1;
        if (niter > INV_JACOBI_SN_MAXITER)
          throw std::invalid_argument("Landen transformation not converging");
    }

    T K = product(1 + ks.slice(1)) * c_pi<T> / 2.;

    univector<complex<T>> wns(1, w);
    univector<T> kn    = ks.slice(0, ks.size()-1);
    univector<T> knext = ks.slice(1, ks.size());

    for (int i = 0, ii = kn.size(); i < ii; i++)
    {
        complex<T> wn    = wns[wns.size()-1];
        complex<T> wnext = (2 * wn / (1 + knext[i]) / (1 + cmpl(kn[i] * wn) /* @WARN double precision limitation in cmpl operation for i = 1 */));
        wns.push_back(wnext); // @TODO : Precision issue to fix.. https://www.sciencedirect.com/science/article/pii/S0010465515001733
    }

    complex<T> u = 2 / c_pi<T> * asin(wns[wns.size()-1]);
    
    return K * u;
}


/*** Real inverse Jacobian sc, with complementary modulus

    References
    ----------
    # noqa: E501
    .. [1] https://functions.wolfram.com/EllipticFunctions/JacobiSC/introductions/JacobiPQs/ShowAll.html,
       "Representations through other Jacobi functions"

    ***/

template <typename T>
T inv_jacobi_sc1(T w, T m)
{
    complex<T> zcomplex = inv_jacobi_sn(complex<T> (0, w), m);

    const T _epsilon = 1e-14;
    if (abs(zcomplex.real()) > _epsilon)
      throw std::invalid_argument("pure imaginary number expected");

    return zcomplex.imag();
}

template <typename T>
T elliptic_k (T m)
{
    //
    // P,Q approximation uses Cephes library 
    // (https://github.com/deepmind/torch-cephes/blob/master/cephes/ellf/ellpk.c)
    double MACHEP =  1.11022302462515654042E-16;   /* 2**-53 */
    // double MACHEP =  1.38777878078144567553E-17;   /* 2**-56 */

    static double P[] =
    {
        1.37982864606273237150E-4,
        2.28025724005875567385E-3,
        7.97404013220415179367E-3,
        9.85821379021226008714E-3,
        6.87489687449949877925E-3,
        6.18901033637687613229E-3,
        8.79078273952743772254E-3,
        1.49380448916805252718E-2,
        3.08851465246711995998E-2,
        9.65735902811690126535E-2,
        1.38629436111989062502E0
    };

    static double Q[] =
    {
        2.94078955048598507511E-5,
        9.14184723865917226571E-4,
        5.94058303753167793257E-3,
        1.54850516649762399335E-2,
        2.39089602715924892727E-2,
        3.01204715227604046988E-2,
        3.73774314173823228969E-2,
        4.88280347570998239232E-2,
        7.03124996963957469739E-2,
        1.24999999999870820058E-1,
        4.99999999999999999821E-1
    };

    static double C1 = 1.3862943611198906188E0; /* log(4) */

    double p = 1.-m;
    if( p < 0 || p > 1.0) return 0;
    if( p > MACHEP ) return cephes_polevl(p, P, 10) - log(p)*cephes_polevl(p, Q, 10);
    if( p == 0 ) return NAN;
    
    return C1 - 0.5 * log(p);
}

template <typename T>
T elliptic_km1 (T p) { return elliptic_k(1.-p); }

template <typename T>
T elliptic_deg ( int N, T m1)
{
    T ellipk   = elliptic_k(m1);
    T ellipkm1 = elliptic_km1(m1);
    
    T q1 = exp(-c_pi<T> / ellipk);
    T q  = pow(q1, ellipkm1/N);

    int const _ELLIPDEG_MMAX = 7;
    static T mnum[_ELLIPDEG_MMAX+1] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    static T mden[_ELLIPDEG_MMAX+1] = { 1, 2, 3, 4, 5, 6, 7, 8 };

    T num = 0;
    T den = 0;
    for (int i = 0; i < _ELLIPDEG_MMAX; i++ )
    {
        num += pow(q, (mnum[i] * (mnum[i]+1)));
        den += pow(q, (mden[i]*mden[i]));
    }

    return 16 * q * pow(num / (1 + 2*den), 4);
}

}
}