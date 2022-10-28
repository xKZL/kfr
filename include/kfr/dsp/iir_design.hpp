/** @addtogroup biquad
 *  @{
 */
/*
  Copyright (C) 2016 D Levin (https://www.kfrlib.com)
  This file is part of KFR

  KFR is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  KFR is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with KFR.

  If GPL is not suitable for your project, you must purchase a commercial license to use KFR.
  Buying a commercial license is mandatory as soon as you develop commercial activities without
  disclosing the source code of your own applications.
  See https://www.kfrlib.com for details.
 */
#pragma once

#include "../base/filter.hpp"
#include "../base/pointer.hpp"
#include "../math/hyperbolic.hpp"
#include "../simd/impl/function.hpp"
#include "../simd/operators.hpp"
#include "../simd/vec.hpp"
#include "../testo/assert.hpp"
#include "biquad_design.hpp"
#include "../math/elliptic.hpp"

#include <cmath>
#include <iostream>

// using namespace std;

namespace kfr
{

template <typename T>
struct zpk
{
    univector<complex<T>> z;
    univector<complex<T>> p;
    T k;
};

inline namespace CMT_ARCH_NAME
{
template <typename T>
KFR_FUNCTION double chebyshev_order(std::vector<T> wp, std::vector<T> ws, identity<T> rp, identity<T> rs)
{
    if(wp.size() != ws.size())
      throw std::invalid_argument("bandpass `wp` and bandstop `ws` coefficient must be of same size");    
    if(wp.size() > 2)
      throw std::invalid_argument("too many coefficients provided for `wp` and `ws`");

    //
    // Deduce filter type based on `wp` and `ws` parameters
    int filter_type = -1;
    if(wp.size() == 1) {

        if(wp[0] < ws[0]) filter_type = 0; // Lowpass
        if(wp[0] > ws[0]) filter_type = 1; // Highpass 

    } else if(wp.size() == 2) {

        if(wp[0] > wp[1]) {
            T swp = wp[0];
            wp[0] = wp[1]; wp[1] = swp;
        }

        if(ws[0] > ws[1]) {
            T swp = ws[0];
            ws[0] = ws[1]; ws[1] = swp;
        }

        if(wp[0] < ws[0] && wp[1] > ws[1]) filter_type = 2; // Bandpass
        if(wp[0] > ws[0] && wp[1] < ws[1]) filter_type = 3; // Bandstop
    }

    wp[0] = tan(c_pi<T> * wp[0] / 2.0);
    wp[1] = tan(c_pi<T> * wp[1] / 2.0);
    ws[0] = tan(c_pi<T> * ws[0] / 2.0);
    ws[1] = tan(c_pi<T> * ws[1] / 2.0);

    if(filter_type < 0)
      throw std::invalid_argument("cannot determine filter type from the coefficient provided");

    //
    // Deduce filter type based on `wp` and `ws` parameters
    std::vector<T> nat;
    if(filter_type == 0) {
    
        nat.push_back(ws[0]/wp[0]);
    
    } else if(filter_type == 1) {

        nat.push_back(wp[0]/ws[0]);

    } else if(filter_type == 2) {

        std::vector<T>  yp0, yp1;
        const T _WEPSILON = 1e-12;

        std::vector<T> _wp = {wp[0], wp[1]};
        
        const int _MAX_STEP = 500;  // # iter limitation
        const T _XTOL       = 1e-4; //     df limitation
        double df0 = max((ws[0] - _WEPSILON - wp[0])/_MAX_STEP, _XTOL);

        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[0] = wp[0] + df0*iStep;

            std::vector<T> nat;
                           nat.push_back( ws[0] * (_wp[0] - _wp[1]) / ( pow(ws[0],2) - _wp[0] * _wp[1] ));
                           nat.push_back( ws[1] * (_wp[0] - _wp[1]) / ( pow(ws[1],2) - _wp[0] * _wp[1] ));

            T _rs  = pow(10, 0.1 * rs);
            T _rp  = pow(10, 0.1 * rp);
            T _nat = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

            yp1.push_back(acosh(sqrt((_rs - 1.0) / (_rp - 1.0))) / acosh(_nat));
        }

        _wp[0] = wp[0] + df0 * (std::min_element(yp0.begin(), yp0.end()) - yp0.begin());

        double df1 = max((wp[1] - (ws[1] + _WEPSILON))/_MAX_STEP, _XTOL);
        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[1] = ws[1] + _WEPSILON + df1*iStep;
            
            std::vector<T> nat;
                           nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
                           nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));

            T _rs  = pow(10, 0.1 * rs);
            T _rp  = pow(10, 0.1 * rp);
            T _nat = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

            yp1.push_back(acosh(sqrt((_rs - 1.0) / (_rp - 1.0))) / acosh(_nat));
        }

        _wp[1] = ws[1] + _WEPSILON + df1 * (std::min_element(yp1.begin(), yp1.end()) - yp1.begin());

        nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
        nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));

    } else if(filter_type == 3) {

        nat.push_back( (pow(ws[0],2) - wp[0] * wp[1]) / (ws[0] * (wp[0] - wp[1])) );
        nat.push_back( (pow(ws[1],2) - wp[0] * wp[1]) / (ws[1] * (wp[0] - wp[1])) );
    }

    T _rs = pow(10, 0.1 * rs);
    T _rp = pow(10, 0.1 * rp);
    T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);
    return ceil(acosh(sqrt((_rs - 1.0) / (_rp - 1.0))) / acosh(_nat));
}

template <typename T>
KFR_FUNCTION zpk<T> chebyshev1(int N, identity<T> rp)
{
    if (N <= 0)
    {
        return { {}, {}, 1 };
    }

    T eps = sqrt(exp10(0.1 * rp) - 1.0);
    T mu  = 1.0 / N * std::asinh(1 / eps);

    univector<T> m = linspace(-N + 1, N + 1, N, false, true);

    univector<T> theta      = c_pi<T> * m / (2 * N);
    univector<complex<T>> p = -csinh(make_complex(mu, theta));
    T k                     = product(-p).real();
    if (N % 2 == 0)
        k = k / sqrt(1.0 + eps * eps);
    return { {}, std::move(p), k };
}

template <typename T>
KFR_FUNCTION zpk<T> chebyshev2(int N, identity<T> rs)
{
    if (N <= 0)
    {
        return { {}, {}, 1 };
    }

    T de = 1.0 / sqrt(exp10(0.1 * rs) - 1);
    T mu = std::asinh(1.0 / de) / N;

    univector<T> m;

    if (N % 2)
    {
        m = concatenate(linspace(-N + 1, -2, N / 2, true, true), linspace(2, N - 1, N / 2, true, true));
    }
    else
    {
        m = linspace(-N + 1, N + 1, N, false, true);
    }

    univector<complex<T>> z = -cconj(complex<T>(0, 1) / sin(m * c_pi<T> / (2.0 * N)));

    univector<complex<T>> p =
        -cexp(complex<T>(0, 1) * c_pi<T> * linspace(-N + 1, N + 1, N, false, true) / (2 * N));
    p = make_complex(sinh(mu) * real(p), cosh(mu) * imag(p));
    p = 1.0 / p;

    T k = (product(-p) / product(-z)).real();

    return { std::move(z), std::move(p), k };
}

template <typename T>
KFR_FUNCTION double elliptic_order(std::vector<T> wp, std::vector<T> ws, identity<T> rp, identity<T> rs)
{
    if(wp.size() != ws.size())
      throw std::invalid_argument("bandpass `wp` and bandstop `ws` coefficient must be of same size");    
    if(wp.size() > 2)
      throw std::invalid_argument("too many coefficients provided for `wp` and `ws`");
    
    //
    // Deduce filter type based on `wp` and `ws` parameters
    int filter_type = -1;
    if(wp.size() == 1) {

        if(wp[0] < ws[0]) filter_type = 0; // Lowpass
        if(wp[0] > ws[0]) filter_type = 1; // Highpass 

    } else if(wp.size() == 2) {

        if(wp[0] > wp[1]) {
            T swp = wp[0];
            wp[0] = wp[1]; wp[1] = swp;
        }

        if(ws[0] > ws[1]) {
            T swp = ws[0];
            ws[0] = ws[1]; ws[1] = swp;
        }

        if(wp[0] < ws[0] && wp[1] > ws[1]) filter_type = 2; // Bandpass
        if(wp[0] > ws[0] && wp[1] < ws[1]) filter_type = 3; // Bandstop
    }

    wp[0] = tan(c_pi<T> * wp[0] / 2.0);
    wp[1] = tan(c_pi<T> * wp[1] / 2.0);
    ws[0] = tan(c_pi<T> * ws[0] / 2.0);
    ws[1] = tan(c_pi<T> * ws[1] / 2.0);

    if(filter_type < 0)
      throw std::invalid_argument("cannot determine filter type from the coefficient provided");

    //
    // Deduce filter type based on `wp` and `ws` parameters
    std::vector<T> nat;
    if(filter_type == 0) {
    
        nat.push_back(ws[0]/wp[0]);
    
    } else if(filter_type == 1) {

        nat.push_back(wp[0]/ws[0]);

    } else if(filter_type == 2) {

        std::vector<T>  yp0, yp1;
        const T _WEPSILON = 1e-12;

        std::vector<T> _wp = {wp[0], wp[1]};
        
        const int _MAX_STEP = 500;  // # iter limitation
        const T _XTOL       = 1e-4; //     df limitation
        double df0 = max((ws[0] - _WEPSILON - wp[0])/_MAX_STEP, _XTOL);

        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[0] = wp[0] + df0*iStep;

            std::vector<T> nat;
                           nat.push_back( ws[0] * (_wp[0] - _wp[1]) / ( pow(ws[0],2) - _wp[0] * _wp[1] ));
                           nat.push_back( ws[1] * (_wp[0] - _wp[1]) / ( pow(ws[1],2) - _wp[0] * _wp[1] ));

            T _rs = pow(10, 0.1 * rs);
            T _rp = pow(10, 0.1 * rp);

            T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);
            T _arg0 = 1. / _nat;
            T _arg1 = sqrt((_rp - 1.0) / (_rs - 1.0));
            
            yp0.push_back(
                (elliptic_k(    pow(_arg0,2)) * elliptic_k(1 - pow(_arg1,2))) / 
                (elliptic_k(1 - pow(_arg0,2)) * elliptic_k(    pow(_arg1,2)))
            );
        }

        _wp[0] = wp[0] + df0 * (std::min_element(yp0.begin(), yp0.end()) - yp0.begin());

        double df1 = max((wp[1] - (ws[1] + _WEPSILON))/_MAX_STEP, _XTOL);
        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[1] = ws[1] + _WEPSILON + df1*iStep;
            
            std::vector<T> nat;
                           nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
                           nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));
            
            T _rs = pow(10, 0.1 * rs);
            T _rp = pow(10, 0.1 * rp);

            T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);
            T _arg0 = 1. / _nat;
            T _arg1 = sqrt((_rp - 1.0) / (_rs - 1.0));

            yp1.push_back(
                (elliptic_k(    pow(_arg0,2)) * elliptic_k(1 - pow(_arg1,2))) / 
                (elliptic_k(1 - pow(_arg0,2)) * elliptic_k(    pow(_arg1,2)))
            );
        }

        _wp[1] = ws[1] + _WEPSILON + df1 * (std::min_element(yp1.begin(), yp1.end()) - yp1.begin());

        nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
        nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));

    } else if(filter_type == 3) {

        nat.push_back( (pow(ws[0],2) - wp[0] * wp[1]) / (ws[0] * (wp[0] - wp[1])) );
        nat.push_back( (pow(ws[1],2) - wp[0] * wp[1]) / (ws[1] * (wp[0] - wp[1])) );
    }

    T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

    T arg0 = 1. / _nat;
    T arg1_sq = (exp10(0.1*rp) - 1) / (exp10(0.1*rs) - 1);

    return ceil(
        (elliptic_k  ( pow(arg0, 2) ) * elliptic_km1(arg1_sq)) / 
        (elliptic_km1( pow(arg0, 2) ) * elliptic_k  (arg1_sq))
    );
}

template <typename T>
KFR_FUNCTION zpk<T> elliptic(int N, identity<T> rp, identity<T> rs)
{
    if (N == 0)
    {
        return { {}, {}, exp10(-rp/20) };
    }

    if (N == 1)
    {
        univector<complex<T>> p = {-csqrt(complex<T>(1, 0) / (exp10(0.1*rp) - 1) )};
        T k = -p[0].real();
        return { {}, p, k };
    }

    double eps_sq = exp10(0.1*rp) - 1;
    double eps = sqrt(eps_sq);

    double ck1_sq = eps_sq / (exp10(0.1*rs) - 1);
    if(ck1_sq == 0)
      throw std::invalid_argument("Cannot design a filter with given rp and rs specifications.");

    double ellipk = elliptic_k(ck1_sq);

    double m      = elliptic_deg(N, ck1_sq);
    double capk   = elliptic_k(m);

    univector<complex<T>> z1, z2;
    univector<T> j;
    for(int i = 1 - N % 2; i <= N; i = i+2)
      j.push_back(i);

    const double _EPSILON = 2e-16;
    int jj = j.size();

    univector<T> s(jj);
    univector<T> c(jj);
    univector<T> d(jj);

    for(int n = 0; n < jj; n++)
    {
        sncndn(j[n] * capk / N, m, s[n], c[n], d[n]);

        if(abs(s[n]) > _EPSILON) 
        {
            z1.push_back( complex<T>(0,  1.0 / (sqrt(m) * s[n])) );
            z2.push_back( complex<T>(0,  1.0 / (sqrt(m) * s[n])) );
        }
    }

    univector<complex<T>> z = concatenate(z1, cconj(z2));

    T sv, cv, dv;
    T r = inv_jacobi_sc1(1. / eps, ck1_sq);
    T v0 = capk * r / (N * ellipk);
    sncndn(v0, 1 - m, sv, cv, dv);

    univector<T> mag2;
    univector<complex<T>> p1, p2;
    for(int n = 0, jj = j.size(); n < jj; n++)
    {
      p1.push_back(- complex<T>(c[n] * d[n] * sv * cv, s[n] * dv) / (1 - pow((d[n] * sv), 2.0)));
      mag2.push_back( p1[n].real()*p1[n].real() + p1[n].imag()*p1[n].imag() );
    }

    for(int n = 0, jj = j.size(); n < jj; n++)
    {
      if(N % 2 == 0 )
        p2.push_back( - complex<T>(c[n] * d[n] * sv * cv, s[n] * dv) / (1 - pow((d[n] * sv), 2.0)) );
      else if(abs(p1[n].imag()) > _EPSILON * sqrt(sum(mag2)))
        p2.push_back( - complex<T>(c[n] * d[n] * sv * cv, s[n] * dv) / (1 - pow((d[n] * sv), 2.0)) );
    }

    univector<complex<T>> p = concatenate(p1, cconj(p2));
    T k = (product(-p) / product(-z)).real();
    
    if (N % 2 == 0)
        k = k / sqrt((1 + eps_sq));
    
    return { std::move(z), std::move(p), k };
}

template <typename T>
KFR_FUNCTION double butterworth_order(std::vector<T> wp, std::vector<T> ws, identity<T> rp, identity<T> rs)
{
    if(wp.size() != ws.size())
      throw std::invalid_argument("bandpass `wp` and bandstop `ws` coefficient must be of same size");    
    if(wp.size() > 2)
      throw std::invalid_argument("too many coefficients provided for `wp` and `ws`");
    
    //
    // Deduce filter type based on `wp` and `ws` parameters
    int filter_type = -1;
    if(wp.size() == 1) {

        if(wp[0] < ws[0]) filter_type = 0; // Lowpass
        if(wp[0] > ws[0]) filter_type = 1; // Highpass 

    } else if(wp.size() == 2) {

        if(wp[0] > wp[1]) {
            T swp = wp[0];
            wp[0] = wp[1]; wp[1] = swp;
        }

        if(ws[0] > ws[1]) {
            T swp = ws[0];
            ws[0] = ws[1]; ws[1] = swp;
        }

        if(wp[0] < ws[0] && wp[1] > ws[1]) filter_type = 2; // Bandpass
        if(wp[0] > ws[0] && wp[1] < ws[1]) filter_type = 3; // Bandstop
    }

    wp[0] = tan(c_pi<T> * wp[0] / 2.0);
    wp[1] = tan(c_pi<T> * wp[1] / 2.0);
    ws[0] = tan(c_pi<T> * ws[0] / 2.0);
    ws[1] = tan(c_pi<T> * ws[1] / 2.0);

    if(filter_type < 0)
      throw std::invalid_argument("cannot determine filter type from the coefficient provided");

    //
    // Deduce filter type based on `wp` and `ws` parameters
    std::vector<T> nat;
    if(filter_type == 0) {
    
        nat.push_back(ws[0]/wp[0]);
    
    } else if(filter_type == 1) {

        nat.push_back(wp[0]/ws[0]);

    } else if(filter_type == 2) {

        std::vector<T>  yp0, yp1;
        const T _WEPSILON = 1e-12;

        std::vector<T> _wp = {wp[0], wp[1]};
        
        const int _MAX_STEP = 500;  // # iter limitation
        const T _XTOL       = 1e-4; //     df limitation
        double df0 = max((ws[0] - _WEPSILON - wp[0])/_MAX_STEP, _XTOL);

        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[0] = wp[0] + df0*iStep;

            std::vector<T> nat;
                      nat.push_back( ws[0] * (_wp[0] - _wp[1]) / ( pow(ws[0],2) - _wp[0] * _wp[1] ));
                      nat.push_back( ws[1] * (_wp[0] - _wp[1]) / ( pow(ws[1],2) - _wp[0] * _wp[1] ));

            T _rs = pow(10, 0.1 * rs);
            T _rp = pow(10, 0.1 * rp);
            T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

            yp0.push_back(log10((_rs - 1.0) / (_rp - 1.0)) / (2 * log10(_nat)));
        }

        _wp[0] = wp[0] + df0 * (std::min_element(yp0.begin(), yp0.end()) - yp0.begin());

        double df1 = max((wp[1] - (ws[1] + _WEPSILON))/_MAX_STEP, _XTOL);
        for(int iStep = 0; iStep <= _MAX_STEP; iStep++) {
        
            _wp[1] = ws[1] + _WEPSILON + df1*iStep;
            
            std::vector<T> nat;
                      nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
                      nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));

            T _rs = pow(10, 0.1 * rs);
            T _rp = pow(10, 0.1 * rp);
            T _nat  = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

            yp1.push_back(log10((_rs - 1.0) / (_rp - 1.0)) / (2 * log10(_nat)));
        }

        _wp[1] = ws[1] + _WEPSILON + df1 * (std::min_element(yp1.begin(), yp1.end()) - yp1.begin());

        nat.push_back(ws[0] * (_wp[0] - _wp[1]) / (pow(ws[0],2) - _wp[0] * _wp[1]));
        nat.push_back(ws[1] * (_wp[0] - _wp[1]) / (pow(ws[1],2) - _wp[0] * _wp[1]));

    } else if(filter_type == 3) {

        nat.push_back( (pow(ws[0],2) - wp[0] * wp[1]) / (ws[0] * (wp[0] - wp[1])) );
        nat.push_back( (pow(ws[1],2) - wp[0] * wp[1]) / (ws[1] * (wp[0] - wp[1])) );
    }

    T _rs  = pow(10, 0.1 * rs);
    T _rp  = pow(10, 0.1 * rp);
    T _nat = nat.size() > 1 ? min(abs(nat[0]), abs(nat[1])) : abs(nat[0]);

    return ceil(log10((_rs - 1.0) / (_rp - 1.0)) / (2 * log10(_nat)));
}

template <typename T>
KFR_FUNCTION zpk<T> butterworth(int N)
{
    switch (N)
    {
    case 1:
        return { {}, { complex<T>(-1., +0.) }, 1 };
    case 2:
        return { {},
                 { complex<T>(-0.7071067811865476, -0.7071067811865476),
                   complex<T>(-0.7071067811865476, +0.7071067811865476) },
                 1 };
    case 3:
        return { {},
                 { complex<T>(-0.5000000000000001, -0.8660254037844386), complex<T>(-1., +0.),
                   complex<T>(-0.5000000000000001, +0.8660254037844386) },
                 1 };
    case 4:
        return { {},
                 { complex<T>(-0.38268343236508984, -0.9238795325112867),
                   complex<T>(-0.9238795325112867, -0.3826834323650898),
                   complex<T>(-0.9238795325112867, +0.3826834323650898),
                   complex<T>(-0.38268343236508984, +0.9238795325112867) },
                 1 };
    case 5:
        return { {},
                 { complex<T>(-0.30901699437494745, -0.9510565162951535),
                   complex<T>(-0.8090169943749475, -0.5877852522924731), complex<T>(-1., +0.),
                   complex<T>(-0.8090169943749475, +0.5877852522924731),
                   complex<T>(-0.30901699437494745, +0.9510565162951535) },
                 1 };
    case 6:
        return { {},
                 { complex<T>(-0.25881904510252096, -0.9659258262890682),
                   complex<T>(-0.7071067811865476, -0.7071067811865476),
                   complex<T>(-0.9659258262890683, -0.25881904510252074),
                   complex<T>(-0.9659258262890683, +0.25881904510252074),
                   complex<T>(-0.7071067811865476, +0.7071067811865476),
                   complex<T>(-0.25881904510252096, +0.9659258262890682) },
                 1 };
    case 7:
        return { {},
                 { complex<T>(-0.22252093395631445, -0.9749279121818236),
                   complex<T>(-0.6234898018587336, -0.7818314824680298),
                   complex<T>(-0.9009688679024191, -0.4338837391175581), complex<T>(-1., +0.),
                   complex<T>(-0.9009688679024191, +0.4338837391175581),
                   complex<T>(-0.6234898018587336, +0.7818314824680298),
                   complex<T>(-0.22252093395631445, +0.9749279121818236) },
                 1 };
    case 8:
        return { {},
                 { complex<T>(-0.19509032201612833, -0.9807852804032304),
                   complex<T>(-0.5555702330196023, -0.8314696123025452),
                   complex<T>(-0.8314696123025452, -0.5555702330196022),
                   complex<T>(-0.9807852804032304, -0.19509032201612825),
                   complex<T>(-0.9807852804032304, +0.19509032201612825),
                   complex<T>(-0.8314696123025452, +0.5555702330196022),
                   complex<T>(-0.5555702330196023, +0.8314696123025452),
                   complex<T>(-0.19509032201612833, +0.9807852804032304) },
                 1 };
    case 9:
        return { {},
                 { complex<T>(-0.17364817766693041, -0.984807753012208),
                   complex<T>(-0.5000000000000001, -0.8660254037844386),
                   complex<T>(-0.766044443118978, -0.6427876096865393),
                   complex<T>(-0.9396926207859084, -0.3420201433256687), complex<T>(-1., +0.),
                   complex<T>(-0.9396926207859084, +0.3420201433256687),
                   complex<T>(-0.766044443118978, +0.6427876096865393),
                   complex<T>(-0.5000000000000001, +0.8660254037844386),
                   complex<T>(-0.17364817766693041, +0.984807753012208) },
                 1 };
    case 10:
        return { {},
                 { complex<T>(-0.15643446504023092, -0.9876883405951378),
                   complex<T>(-0.4539904997395468, -0.8910065241883678),
                   complex<T>(-0.7071067811865476, -0.7071067811865476),
                   complex<T>(-0.8910065241883679, -0.45399049973954675),
                   complex<T>(-0.9876883405951378, -0.15643446504023087),
                   complex<T>(-0.9876883405951378, +0.15643446504023087),
                   complex<T>(-0.8910065241883679, +0.45399049973954675),
                   complex<T>(-0.7071067811865476, +0.7071067811865476),
                   complex<T>(-0.4539904997395468, +0.8910065241883678),
                   complex<T>(-0.15643446504023092, +0.9876883405951378) },
                 1 };
    case 11:
        return { {},
                 { complex<T>(-0.14231483827328512, -0.9898214418809327),
                   complex<T>(-0.41541501300188644, -0.9096319953545183),
                   complex<T>(-0.654860733945285, -0.7557495743542583),
                   complex<T>(-0.8412535328311812, -0.5406408174555976),
                   complex<T>(-0.9594929736144974, -0.28173255684142967), complex<T>(-1., +0.),
                   complex<T>(-0.9594929736144974, +0.28173255684142967),
                   complex<T>(-0.8412535328311812, +0.5406408174555976),
                   complex<T>(-0.654860733945285, +0.7557495743542583),
                   complex<T>(-0.41541501300188644, +0.9096319953545183),
                   complex<T>(-0.14231483827328512, +0.9898214418809327) },
                 1 };
    case 12:
        return { {},
                 { complex<T>(-0.13052619222005193, -0.9914448613738104),
                   complex<T>(-0.38268343236508984, -0.9238795325112867),
                   complex<T>(-0.6087614290087207, -0.7933533402912352),
                   complex<T>(-0.7933533402912353, -0.6087614290087205),
                   complex<T>(-0.9238795325112867, -0.3826834323650898),
                   complex<T>(-0.9914448613738104, -0.13052619222005157),
                   complex<T>(-0.9914448613738104, +0.13052619222005157),
                   complex<T>(-0.9238795325112867, +0.3826834323650898),
                   complex<T>(-0.7933533402912353, +0.6087614290087205),
                   complex<T>(-0.6087614290087207, +0.7933533402912352),
                   complex<T>(-0.38268343236508984, +0.9238795325112867),
                   complex<T>(-0.13052619222005193, +0.9914448613738104) },
                 1 };
    case 13:
        return { {},
                 { complex<T>(-0.120536680255323, -0.992708874098054),
                   complex<T>(-0.35460488704253557, -0.9350162426854148),
                   complex<T>(-0.5680647467311558, -0.8229838658936565),
                   complex<T>(-0.7485107481711011, -0.6631226582407952),
                   complex<T>(-0.8854560256532099, -0.46472317204376856),
                   complex<T>(-0.970941817426052, -0.23931566428755777), complex<T>(-1., +0.),
                   complex<T>(-0.970941817426052, +0.23931566428755777),
                   complex<T>(-0.8854560256532099, +0.46472317204376856),
                   complex<T>(-0.7485107481711011, +0.6631226582407952),
                   complex<T>(-0.5680647467311558, +0.8229838658936565),
                   complex<T>(-0.35460488704253557, +0.9350162426854148),
                   complex<T>(-0.120536680255323, +0.992708874098054) },
                 1 };
    case 14:
        return { {},
                 { complex<T>(-0.11196447610330791, -0.9937122098932426),
                   complex<T>(-0.3302790619551673, -0.9438833303083675),
                   complex<T>(-0.5320320765153366, -0.8467241992282841),
                   complex<T>(-0.7071067811865476, -0.7071067811865476),
                   complex<T>(-0.8467241992282842, -0.5320320765153366),
                   complex<T>(-0.9438833303083676, -0.3302790619551671),
                   complex<T>(-0.9937122098932426, -0.11196447610330786),
                   complex<T>(-0.9937122098932426, +0.11196447610330786),
                   complex<T>(-0.9438833303083676, +0.3302790619551671),
                   complex<T>(-0.8467241992282842, +0.5320320765153366),
                   complex<T>(-0.7071067811865476, +0.7071067811865476),
                   complex<T>(-0.5320320765153366, +0.8467241992282841),
                   complex<T>(-0.3302790619551673, +0.9438833303083675),
                   complex<T>(-0.11196447610330791, +0.9937122098932426) },
                 1 };
    case 15:
        return { {},
                 { complex<T>(-0.10452846326765346, -0.9945218953682733),
                   complex<T>(-0.30901699437494745, -0.9510565162951535),
                   complex<T>(-0.5000000000000001, -0.8660254037844386),
                   complex<T>(-0.6691306063588582, -0.7431448254773941),
                   complex<T>(-0.8090169943749475, -0.5877852522924731),
                   complex<T>(-0.9135454576426009, -0.40673664307580015),
                   complex<T>(-0.9781476007338057, -0.20791169081775931), complex<T>(-1., +0.),
                   complex<T>(-0.9781476007338057, +0.20791169081775931),
                   complex<T>(-0.9135454576426009, +0.40673664307580015),
                   complex<T>(-0.8090169943749475, +0.5877852522924731),
                   complex<T>(-0.6691306063588582, +0.7431448254773941),
                   complex<T>(-0.5000000000000001, +0.8660254037844386),
                   complex<T>(-0.30901699437494745, +0.9510565162951535),
                   complex<T>(-0.10452846326765346, +0.9945218953682733) },
                 1 };
    case 16:
        return { {},
                 { complex<T>(-0.09801714032956077, -0.9951847266721968),
                   complex<T>(-0.29028467725446233, -0.9569403357322089),
                   complex<T>(-0.4713967368259978, -0.8819212643483549),
                   complex<T>(-0.6343932841636455, -0.7730104533627369),
                   complex<T>(-0.773010453362737, -0.6343932841636455),
                   complex<T>(-0.881921264348355, -0.47139673682599764),
                   complex<T>(-0.9569403357322088, -0.29028467725446233),
                   complex<T>(-0.9951847266721969, -0.0980171403295606),
                   complex<T>(-0.9951847266721969, +0.0980171403295606),
                   complex<T>(-0.9569403357322088, +0.29028467725446233),
                   complex<T>(-0.881921264348355, +0.47139673682599764),
                   complex<T>(-0.773010453362737, +0.6343932841636455),
                   complex<T>(-0.6343932841636455, +0.7730104533627369),
                   complex<T>(-0.4713967368259978, +0.8819212643483549),
                   complex<T>(-0.29028467725446233, +0.9569403357322089),
                   complex<T>(-0.09801714032956077, +0.9951847266721968) },
                 1 };
    case 17:
        return { {},
                 { complex<T>(-0.09226835946330202, -0.9957341762950345),
                   complex<T>(-0.273662990072083, -0.961825643172819),
                   complex<T>(-0.4457383557765383, -0.8951632913550623),
                   complex<T>(-0.6026346363792565, -0.7980172272802395),
                   complex<T>(-0.7390089172206591, -0.6736956436465572),
                   complex<T>(-0.8502171357296142, -0.5264321628773557),
                   complex<T>(-0.9324722294043558, -0.3612416661871529),
                   complex<T>(-0.9829730996839018, -0.18374951781657034), complex<T>(-1., +0.),
                   complex<T>(-0.9829730996839018, +0.18374951781657034),
                   complex<T>(-0.9324722294043558, +0.3612416661871529),
                   complex<T>(-0.8502171357296142, +0.5264321628773557),
                   complex<T>(-0.7390089172206591, +0.6736956436465572),
                   complex<T>(-0.6026346363792565, +0.7980172272802395),
                   complex<T>(-0.4457383557765383, +0.8951632913550623),
                   complex<T>(-0.273662990072083, +0.961825643172819),
                   complex<T>(-0.09226835946330202, +0.9957341762950345) },
                 1 };
    case 18:
        return { {},
                 { complex<T>(-0.08715574274765814, -0.9961946980917455),
                   complex<T>(-0.25881904510252096, -0.9659258262890682),
                   complex<T>(-0.42261826174069944, -0.9063077870366499),
                   complex<T>(-0.5735764363510463, -0.8191520442889917),
                   complex<T>(-0.7071067811865476, -0.7071067811865476),
                   complex<T>(-0.8191520442889918, -0.573576436351046),
                   complex<T>(-0.9063077870366499, -0.4226182617406994),
                   complex<T>(-0.9659258262890683, -0.25881904510252074),
                   complex<T>(-0.9961946980917455, -0.08715574274765817),
                   complex<T>(-0.9961946980917455, +0.08715574274765817),
                   complex<T>(-0.9659258262890683, +0.25881904510252074),
                   complex<T>(-0.9063077870366499, +0.4226182617406994),
                   complex<T>(-0.8191520442889918, +0.573576436351046),
                   complex<T>(-0.7071067811865476, +0.7071067811865476),
                   complex<T>(-0.5735764363510463, +0.8191520442889917),
                   complex<T>(-0.42261826174069944, +0.9063077870366499),
                   complex<T>(-0.25881904510252096, +0.9659258262890682),
                   complex<T>(-0.08715574274765814, +0.9961946980917455) },
                 1 };
    case 19:
        return { {},
                 { complex<T>(-0.0825793454723324, -0.9965844930066698),
                   complex<T>(-0.24548548714079924, -0.9694002659393304),
                   complex<T>(-0.40169542465296953, -0.9157733266550574),
                   complex<T>(-0.5469481581224269, -0.8371664782625285),
                   complex<T>(-0.6772815716257411, -0.7357239106731316),
                   complex<T>(-0.7891405093963936, -0.6142127126896678),
                   complex<T>(-0.8794737512064891, -0.4759473930370735),
                   complex<T>(-0.9458172417006346, -0.32469946920468346),
                   complex<T>(-0.9863613034027223, -0.1645945902807339), complex<T>(-1., +0.),
                   complex<T>(-0.9863613034027223, +0.1645945902807339),
                   complex<T>(-0.9458172417006346, +0.32469946920468346),
                   complex<T>(-0.8794737512064891, +0.4759473930370735),
                   complex<T>(-0.7891405093963936, +0.6142127126896678),
                   complex<T>(-0.6772815716257411, +0.7357239106731316),
                   complex<T>(-0.5469481581224269, +0.8371664782625285),
                   complex<T>(-0.40169542465296953, +0.9157733266550574),
                   complex<T>(-0.24548548714079924, +0.9694002659393304),
                   complex<T>(-0.0825793454723324, +0.9965844930066698) },
                 1 };
    case 20:
        return { {},
                 { complex<T>(-0.078459095727845, -0.996917333733128),
                   complex<T>(-0.23344536385590525, -0.9723699203976767),
                   complex<T>(-0.38268343236508984, -0.9238795325112867),
                   complex<T>(-0.5224985647159487, -0.8526401643540923),
                   complex<T>(-0.6494480483301837, -0.7604059656000309),
                   complex<T>(-0.7604059656000309, -0.6494480483301837),
                   complex<T>(-0.8526401643540923, -0.5224985647159488),
                   complex<T>(-0.9238795325112867, -0.3826834323650898),
                   complex<T>(-0.9723699203976766, -0.2334453638559054),
                   complex<T>(-0.996917333733128, -0.07845909572784494),
                   complex<T>(-0.996917333733128, +0.07845909572784494),
                   complex<T>(-0.9723699203976766, +0.2334453638559054),
                   complex<T>(-0.9238795325112867, +0.3826834323650898),
                   complex<T>(-0.8526401643540923, +0.5224985647159488),
                   complex<T>(-0.7604059656000309, +0.6494480483301837),
                   complex<T>(-0.6494480483301837, +0.7604059656000309),
                   complex<T>(-0.5224985647159487, +0.8526401643540923),
                   complex<T>(-0.38268343236508984, +0.9238795325112867),
                   complex<T>(-0.23344536385590525, +0.9723699203976767),
                   complex<T>(-0.078459095727845, +0.996917333733128) },
                 1 };
    case 21:
        return { {},
                 { complex<T>(-0.07473009358642439, -0.9972037971811801),
                   complex<T>(-0.22252093395631445, -0.9749279121818236),
                   complex<T>(-0.3653410243663952, -0.9308737486442041),
                   complex<T>(-0.5000000000000001, -0.8660254037844386),
                   complex<T>(-0.6234898018587336, -0.7818314824680298),
                   complex<T>(-0.7330518718298263, -0.6801727377709194),
                   complex<T>(-0.8262387743159949, -0.563320058063622),
                   complex<T>(-0.9009688679024191, -0.4338837391175581),
                   complex<T>(-0.9555728057861408, -0.29475517441090415),
                   complex<T>(-0.9888308262251285, -0.14904226617617441),
                   complex<T>(-1., +0.),
                   complex<T>(-0.9888308262251285, +0.14904226617617441),
                   complex<T>(-0.9555728057861408, +0.29475517441090415),
                   complex<T>(-0.9009688679024191, +0.4338837391175581),
                   complex<T>(-0.8262387743159949, +0.563320058063622),
                   complex<T>(-0.7330518718298263, +0.6801727377709194),
                   complex<T>(-0.6234898018587336, +0.7818314824680298),
                   complex<T>(-0.5000000000000001, +0.8660254037844386),
                   complex<T>(-0.3653410243663952, +0.9308737486442041),
                   complex<T>(-0.22252093395631445, +0.9749279121818236),
                   complex<T>(-0.07473009358642439, +0.9972037971811801) },
                 1 };
    case 22:
        return { {},
                 { complex<T>(-0.07133918319923235, -0.9974521146102535),
                   complex<T>(-0.21256528955297682, -0.9771468659711595),
                   complex<T>(-0.3494641795990982, -0.9369497249997618),
                   complex<T>(-0.479248986720057, -0.8776789895672555),
                   complex<T>(-0.5992776665113468, -0.8005412409243605),
                   complex<T>(-0.7071067811865477, -0.7071067811865475),
                   complex<T>(-0.8005412409243604, -0.5992776665113468),
                   complex<T>(-0.8776789895672557, -0.4792489867200568),
                   complex<T>(-0.9369497249997617, -0.34946417959909837),
                   complex<T>(-0.9771468659711595, -0.21256528955297668),
                   complex<T>(-0.9974521146102535, -0.07133918319923234),
                   complex<T>(-0.9974521146102535, +0.07133918319923234),
                   complex<T>(-0.9771468659711595, +0.21256528955297668),
                   complex<T>(-0.9369497249997617, +0.34946417959909837),
                   complex<T>(-0.8776789895672557, +0.4792489867200568),
                   complex<T>(-0.8005412409243604, +0.5992776665113468),
                   complex<T>(-0.7071067811865477, +0.7071067811865475),
                   complex<T>(-0.5992776665113468, +0.8005412409243605),
                   complex<T>(-0.479248986720057, +0.8776789895672555),
                   complex<T>(-0.3494641795990982, +0.9369497249997618),
                   complex<T>(-0.21256528955297682, +0.9771468659711595),
                   complex<T>(-0.07133918319923235, +0.9974521146102535) },
                 1 };
    case 23:
        return { {},
                 { complex<T>(-0.06824241336467123, -0.9976687691905392),
                   complex<T>(-0.20345601305263397, -0.9790840876823228),
                   complex<T>(-0.3348796121709863, -0.9422609221188204),
                   complex<T>(-0.4600650377311522, -0.8878852184023752),
                   complex<T>(-0.5766803221148672, -0.816969893010442),
                   complex<T>(-0.6825531432186541, -0.730835964278124),
                   complex<T>(-0.7757112907044199, -0.6310879443260528),
                   complex<T>(-0.8544194045464886, -0.5195839500354336),
                   complex<T>(-0.917211301505453, -0.39840108984624145),
                   complex<T>(-0.9629172873477994, -0.2697967711570243),
                   complex<T>(-0.9906859460363308, -0.1361666490962466),
                   complex<T>(-1., +0.),
                   complex<T>(-0.9906859460363308, +0.1361666490962466),
                   complex<T>(-0.9629172873477994, +0.2697967711570243),
                   complex<T>(-0.917211301505453, +0.39840108984624145),
                   complex<T>(-0.8544194045464886, +0.5195839500354336),
                   complex<T>(-0.7757112907044199, +0.6310879443260528),
                   complex<T>(-0.6825531432186541, +0.730835964278124),
                   complex<T>(-0.5766803221148672, +0.816969893010442),
                   complex<T>(-0.4600650377311522, +0.8878852184023752),
                   complex<T>(-0.3348796121709863, +0.9422609221188204),
                   complex<T>(-0.20345601305263397, +0.9790840876823228),
                   complex<T>(-0.06824241336467123, +0.9976687691905392) },
                 1 };
    case 24:
        return { {},
                 { complex<T>(-0.06540312923014327, -0.9978589232386035),
                   complex<T>(-0.19509032201612833, -0.9807852804032304),
                   complex<T>(-0.3214394653031617, -0.9469301294951056),
                   complex<T>(-0.44228869021900125, -0.8968727415326884),
                   complex<T>(-0.5555702330196024, -0.8314696123025451),
                   complex<T>(-0.6593458151000688, -0.7518398074789774),
                   complex<T>(-0.7518398074789775, -0.6593458151000687),
                   complex<T>(-0.8314696123025452, -0.5555702330196022),
                   complex<T>(-0.8968727415326884, -0.44228869021900125),
                   complex<T>(-0.9469301294951057, -0.32143946530316153),
                   complex<T>(-0.9807852804032304, -0.19509032201612825),
                   complex<T>(-0.9978589232386035, -0.06540312923014306),
                   complex<T>(-0.9978589232386035, +0.06540312923014306),
                   complex<T>(-0.9807852804032304, +0.19509032201612825),
                   complex<T>(-0.9469301294951057, +0.32143946530316153),
                   complex<T>(-0.8968727415326884, +0.44228869021900125),
                   complex<T>(-0.8314696123025452, +0.5555702330196022),
                   complex<T>(-0.7518398074789775, +0.6593458151000687),
                   complex<T>(-0.6593458151000688, +0.7518398074789774),
                   complex<T>(-0.5555702330196024, +0.8314696123025451),
                   complex<T>(-0.44228869021900125, +0.8968727415326884),
                   complex<T>(-0.3214394653031617, +0.9469301294951056),
                   complex<T>(-0.19509032201612833, +0.9807852804032304),
                   complex<T>(-0.06540312923014327, +0.9978589232386035) },
                 1 };
    default:
        return { {}, {}, 1.0 };
    }
}

template <typename T>
KFR_FUNCTION zpk<T> bessel(int N)
{
    switch (N)
    {
    case 0:
        return { {}, {}, 1.0 };
    case 1:
        return { {}, { complex<T>(-1., +0.) }, 1.0 };
    case 2:
        return { {},
                 { complex<T>(-0.8660254037844385, -0.4999999999999999),
                   complex<T>(-0.8660254037844385, +0.4999999999999999) },
                 1.0 };
    case 3:
        return { {},
                 { complex<T>(-0.7456403858480766, -0.7113666249728353), complex<T>(-0.9416000265332067, +0.),
                   complex<T>(-0.7456403858480766, +0.7113666249728353) },
                 1.0 };
    case 4:
        return { {},
                 { complex<T>(-0.6572111716718827, -0.830161435004873),
                   complex<T>(-0.904758796788245, -0.27091873300387465),
                   complex<T>(-0.904758796788245, +0.27091873300387465),
                   complex<T>(-0.6572111716718827, +0.830161435004873) },
                 1.0 };
    case 5:
        return { {},
                 { complex<T>(-0.5905759446119192, -0.9072067564574549),
                   complex<T>(-0.8515536193688396, -0.44271746394433276), complex<T>(-0.92644207738776, +0.),
                   complex<T>(-0.8515536193688396, +0.44271746394433276),
                   complex<T>(-0.5905759446119192, +0.9072067564574549) },
                 1.0 };
    case 6:
        return { {},
                 { complex<T>(-0.5385526816693108, -0.9616876881954276),
                   complex<T>(-0.7996541858328287, -0.5621717346937316),
                   complex<T>(-0.9093906830472273, -0.1856964396793047),
                   complex<T>(-0.9093906830472273, +0.1856964396793047),
                   complex<T>(-0.7996541858328287, +0.5621717346937316),
                   complex<T>(-0.5385526816693108, +0.9616876881954276) },
                 1.0 };
    case 7:
        return { {},
                 { complex<T>(-0.4966917256672317, -1.0025085084544205),
                   complex<T>(-0.7527355434093214, -0.6504696305522553),
                   complex<T>(-0.8800029341523379, -0.32166527623077407), complex<T>(-0.919487155649029, +0.),
                   complex<T>(-0.8800029341523379, +0.32166527623077407),
                   complex<T>(-0.7527355434093214, +0.6504696305522553),
                   complex<T>(-0.4966917256672317, +1.0025085084544205) },
                 1.0 };
    case 8:
        return { {},
                 { complex<T>(-0.46217404125321254, -1.0343886811269012),
                   complex<T>(-0.7111381808485397, -0.71865173141084),
                   complex<T>(-0.8473250802359334, -0.42590175382729345),
                   complex<T>(-0.9096831546652911, -0.1412437976671423),
                   complex<T>(-0.9096831546652911, +0.1412437976671423),
                   complex<T>(-0.8473250802359334, +0.42590175382729345),
                   complex<T>(-0.7111381808485397, +0.71865173141084),
                   complex<T>(-0.46217404125321254, +1.0343886811269012) },
                 1.0 };
    case 9:
        return { {},
                 { complex<T>(-0.43314155615536226, -1.0600736701359301),
                   complex<T>(-0.6743622686854763, -0.7730546212691185),
                   complex<T>(-0.8148021112269014, -0.50858156896315),
                   complex<T>(-0.8911217017079759, -0.25265809345821644),
                   complex<T>(-0.9154957797499037, +0.),
                   complex<T>(-0.8911217017079759, +0.25265809345821644),
                   complex<T>(-0.8148021112269014, +0.50858156896315),
                   complex<T>(-0.6743622686854763, +0.7730546212691185),
                   complex<T>(-0.43314155615536226, +1.0600736701359301) },
                 1.0 };
    case 10:
        return { {},
                 { complex<T>(-0.4083220732868863, -1.0812748428191246),
                   complex<T>(-0.6417513866988322, -0.8175836167191022),
                   complex<T>(-0.7837694413101445, -0.5759147538499948),
                   complex<T>(-0.8688459641284763, -0.34300082337663096),
                   complex<T>(-0.9091347320900505, -0.11395831373355114),
                   complex<T>(-0.9091347320900505, +0.11395831373355114),
                   complex<T>(-0.8688459641284763, +0.34300082337663096),
                   complex<T>(-0.7837694413101445, +0.5759147538499948),
                   complex<T>(-0.6417513866988322, +0.8175836167191022),
                   complex<T>(-0.4083220732868863, +1.0812748428191246) },
                 1.0 };
    case 11:
        return { {},
                 { complex<T>(-0.3868149510055096, -1.0991174667631216),
                   complex<T>(-0.6126871554915198, -0.8547813893314768),
                   complex<T>(-0.7546938934722305, -0.6319150050721849),
                   complex<T>(-0.8453044014712964, -0.4178696917801249),
                   complex<T>(-0.8963656705721169, -0.20804803750710327),
                   complex<T>(-0.9129067244518985, +0.),
                   complex<T>(-0.8963656705721169, +0.20804803750710327),
                   complex<T>(-0.8453044014712964, +0.4178696917801249),
                   complex<T>(-0.7546938934722305, +0.6319150050721849),
                   complex<T>(-0.6126871554915198, +0.8547813893314768),
                   complex<T>(-0.3868149510055096, +1.0991174667631216) },
                 1.0 };
    case 12:
        return { {},
                 { complex<T>(-0.3679640085526314, -1.1143735756415463),
                   complex<T>(-0.5866369321861475, -0.8863772751320727),
                   complex<T>(-0.7276681615395161, -0.6792961178764695),
                   complex<T>(-0.8217296939939074, -0.48102121151006755),
                   complex<T>(-0.880253434201683, -0.2871779503524227),
                   complex<T>(-0.9084478234140686, -0.09550636521345042),
                   complex<T>(-0.9084478234140686, +0.09550636521345042),
                   complex<T>(-0.880253434201683, +0.2871779503524227),
                   complex<T>(-0.8217296939939074, +0.48102121151006755),
                   complex<T>(-0.7276681615395161, +0.6792961178764695),
                   complex<T>(-0.5866369321861475, +0.8863772751320727),
                   complex<T>(-0.3679640085526314, +1.1143735756415463) },
                 1.0 };
    case 13:
        return { {},
                 { complex<T>(-0.35127923233898156, -1.1275915483177048),
                   complex<T>(-0.5631559842430193, -0.9135900338325103),
                   complex<T>(-0.7026234675721276, -0.7199611890171304),
                   complex<T>(-0.7987460692470972, -0.5350752120696802),
                   complex<T>(-0.8625094198260553, -0.35474137311729914),
                   complex<T>(-0.8991314665475196, -0.17683429561610436),
                   complex<T>(-0.9110914665984182, +0.),
                   complex<T>(-0.8991314665475196, +0.17683429561610436),
                   complex<T>(-0.8625094198260553, +0.35474137311729914),
                   complex<T>(-0.7987460692470972, +0.5350752120696802),
                   complex<T>(-0.7026234675721276, +0.7199611890171304),
                   complex<T>(-0.5631559842430193, +0.9135900338325103),
                   complex<T>(-0.35127923233898156, +1.1275915483177048) },
                 1.0 };
    case 14:
        return { {},
                 { complex<T>(-0.3363868224902031, -1.1391722978398593),
                   complex<T>(-0.5418766775112306, -0.9373043683516926),
                   complex<T>(-0.6794256425119227, -0.7552857305042031),
                   complex<T>(-0.7766591387063624, -0.581917067737761),
                   complex<T>(-0.8441199160909851, -0.41316538251026935),
                   complex<T>(-0.8869506674916446, -0.24700791787653334),
                   complex<T>(-0.9077932138396491, -0.08219639941940154),
                   complex<T>(-0.9077932138396491, +0.08219639941940154),
                   complex<T>(-0.8869506674916446, +0.24700791787653334),
                   complex<T>(-0.8441199160909851, +0.41316538251026935),
                   complex<T>(-0.7766591387063624, +0.581917067737761),
                   complex<T>(-0.6794256425119227, +0.7552857305042031),
                   complex<T>(-0.5418766775112306, +0.9373043683516926),
                   complex<T>(-0.3363868224902031, +1.1391722978398593) },
                 1.0 };
    case 15:
        return { {},
                 { complex<T>(-0.3229963059766445, -1.14941615458363),
                   complex<T>(-0.5224954069658334, -0.9581787261092526),
                   complex<T>(-0.6579196593111002, -0.7862895503722519),
                   complex<T>(-0.7556027168970723, -0.6229396358758262),
                   complex<T>(-0.8256631452587148, -0.46423487527343266),
                   complex<T>(-0.8731264620834984, -0.30823524705642674),
                   complex<T>(-0.9006981694176978, -0.1537681197278439), complex<T>(-0.9097482363849062, +0.),
                   complex<T>(-0.9006981694176978, +0.1537681197278439),
                   complex<T>(-0.8731264620834984, +0.30823524705642674),
                   complex<T>(-0.8256631452587148, +0.46423487527343266),
                   complex<T>(-0.7556027168970723, +0.6229396358758262),
                   complex<T>(-0.6579196593111002, +0.7862895503722519),
                   complex<T>(-0.5224954069658334, +0.9581787261092526),
                   complex<T>(-0.3229963059766445, +1.14941615458363) },
                 1.0 };
    case 16:
        return { {},
                 { complex<T>(-0.31087827556453856, -1.1585528411993304),
                   complex<T>(-0.504760644442476, -0.9767137477799086),
                   complex<T>(-0.6379502514039066, -0.8137453537108762),
                   complex<T>(-0.7356166304713119, -0.6591950877860395),
                   complex<T>(-0.8074790293236005, -0.5092933751171799),
                   complex<T>(-0.8584264231521322, -0.3621697271802063),
                   complex<T>(-0.8911723070323642, -0.2167089659900575),
                   complex<T>(-0.9072099595087003, -0.07214211304111734),
                   complex<T>(-0.9072099595087003, +0.07214211304111734),
                   complex<T>(-0.8911723070323642, +0.2167089659900575),
                   complex<T>(-0.8584264231521322, +0.3621697271802063),
                   complex<T>(-0.8074790293236005, +0.5092933751171799),
                   complex<T>(-0.7356166304713119, +0.6591950877860395),
                   complex<T>(-0.6379502514039066, +0.8137453537108762),
                   complex<T>(-0.504760644442476, +0.9767137477799086),
                   complex<T>(-0.31087827556453856, +1.1585528411993304) },
                 1.0 };
    case 17:
        return { {},
                 { complex<T>(-0.29984894599900724, -1.1667612729256673),
                   complex<T>(-0.48846293376727057, -0.9932971956316782),
                   complex<T>(-0.6193710717342136, -0.8382497252826987),
                   complex<T>(-0.7166893842372348, -0.6914936286393606),
                   complex<T>(-0.7897644147799701, -0.5493724405281085),
                   complex<T>(-0.8433414495836128, -0.41007592829100215),
                   complex<T>(-0.8801100704438625, -0.2725347156478803),
                   complex<T>(-0.9016273850787279, -0.13602679951730237),
                   complex<T>(-0.9087141161336388, +0.),
                   complex<T>(-0.9016273850787279, +0.13602679951730237),
                   complex<T>(-0.8801100704438625, +0.2725347156478803),
                   complex<T>(-0.8433414495836128, +0.41007592829100215),
                   complex<T>(-0.7897644147799701, +0.5493724405281085),
                   complex<T>(-0.7166893842372348, +0.6914936286393606),
                   complex<T>(-0.6193710717342136, +0.8382497252826987),
                   complex<T>(-0.48846293376727057, +0.9932971956316782),
                   complex<T>(-0.29984894599900724, +1.1667612729256673) },
                 1.0 };
    case 18:
        return { {},
                 { complex<T>(-0.28975920298804847, -1.1741830106000584),
                   complex<T>(-0.4734268069916154, -1.0082343003148009),
                   complex<T>(-0.6020482668090646, -0.8602708961893666),
                   complex<T>(-0.698782144500527, -0.7204696509726628),
                   complex<T>(-0.7726285030739557, -0.5852778162086639),
                   complex<T>(-0.8281885016242831, -0.45293856978159136),
                   complex<T>(-0.8681095503628832, -0.32242049251632576),
                   complex<T>(-0.8939764278132456, -0.19303746408947586),
                   complex<T>(-0.9067004324162776, -0.06427924106393067),
                   complex<T>(-0.9067004324162776, +0.06427924106393067),
                   complex<T>(-0.8939764278132456, +0.19303746408947586),
                   complex<T>(-0.8681095503628832, +0.32242049251632576),
                   complex<T>(-0.8281885016242831, +0.45293856978159136),
                   complex<T>(-0.7726285030739557, +0.5852778162086639),
                   complex<T>(-0.698782144500527, +0.7204696509726628),
                   complex<T>(-0.6020482668090646, +0.8602708961893666),
                   complex<T>(-0.4734268069916154, +1.0082343003148009),
                   complex<T>(-0.28975920298804847, +1.1741830106000584) },
                 1.0 };
    case 19:
        return { {},
                 { complex<T>(-0.2804866851439361, -1.1809316284532905),
                   complex<T>(-0.4595043449730983, -1.0217687769126707),
                   complex<T>(-0.5858613321217832, -0.8801817131014564),
                   complex<T>(-0.6818424412912442, -0.7466272357947761),
                   complex<T>(-0.7561260971541627, -0.6176483917970176),
                   complex<T>(-0.8131725551578203, -0.491536503556246),
                   complex<T>(-0.8555768765618422, -0.3672925896399872),
                   complex<T>(-0.8849290585034385, -0.24425907575498182),
                   complex<T>(-0.9021937639390656, -0.12195683818720263),
                   complex<T>(-0.9078934217899399, +0.),
                   complex<T>(-0.9021937639390656, +0.12195683818720263),
                   complex<T>(-0.8849290585034385, +0.24425907575498182),
                   complex<T>(-0.8555768765618422, +0.3672925896399872),
                   complex<T>(-0.8131725551578203, +0.491536503556246),
                   complex<T>(-0.7561260971541627, +0.6176483917970176),
                   complex<T>(-0.6818424412912442, +0.7466272357947761),
                   complex<T>(-0.5858613321217832, +0.8801817131014564),
                   complex<T>(-0.4595043449730983, +1.0217687769126707),
                   complex<T>(-0.2804866851439361, +1.1809316284532905) },
                 1.0 };
    case 20:
        return { {},
                 { complex<T>(-0.27192995802516506, -1.187099379810886),
                   complex<T>(-0.44657006982051484, -1.0340977025608422),
                   complex<T>(-0.5707026806915716, -0.8982829066468254),
                   complex<T>(-0.6658120544829932, -0.7703721701100759),
                   complex<T>(-0.7402780309646764, -0.6469975237605227),
                   complex<T>(-0.7984251191290602, -0.526494238881713),
                   complex<T>(-0.8427907479956664, -0.4078917326291931),
                   complex<T>(-0.8749560316673335, -0.2905559296567909),
                   complex<T>(-0.8959150941925766, -0.17403171759187044),
                   complex<T>(-0.9062570115576768, -0.0579617802778495),
                   complex<T>(-0.9062570115576768, +0.0579617802778495),
                   complex<T>(-0.8959150941925766, +0.17403171759187044),
                   complex<T>(-0.8749560316673335, +0.2905559296567909),
                   complex<T>(-0.8427907479956664, +0.4078917326291931),
                   complex<T>(-0.7984251191290602, +0.526494238881713),
                   complex<T>(-0.7402780309646764, +0.6469975237605227),
                   complex<T>(-0.6658120544829932, +0.7703721701100759),
                   complex<T>(-0.5707026806915716, +0.8982829066468254),
                   complex<T>(-0.44657006982051484, +1.0340977025608422),
                   complex<T>(-0.27192995802516506, +1.187099379810886) },
                 1.0 };
    case 21:
        return { {},
                 { complex<T>(-0.2640041595834027, -1.192762031948052),
                   complex<T>(-0.4345168906815268, -1.045382255856986),
                   complex<T>(-0.5564766488918566, -0.9148198405846728),
                   complex<T>(-0.6506315378609466, -0.7920349342629495),
                   complex<T>(-0.7250839687106612, -0.6737426063024383),
                   complex<T>(-0.7840287980408347, -0.5583186348022857),
                   complex<T>(-0.8299435470674444, -0.44481777394079575),
                   complex<T>(-0.8643915813643203, -0.33262585125221866),
                   complex<T>(-0.8883808106664449, -0.221306921508435),
                   complex<T>(-0.9025428073192694, -0.11052525727898564),
                   complex<T>(-0.9072262653142963, +0.),
                   complex<T>(-0.9025428073192694, +0.11052525727898564),
                   complex<T>(-0.8883808106664449, +0.221306921508435),
                   complex<T>(-0.8643915813643203, +0.33262585125221866),
                   complex<T>(-0.8299435470674444, +0.44481777394079575),
                   complex<T>(-0.7840287980408347, +0.5583186348022857),
                   complex<T>(-0.7250839687106612, +0.6737426063024383),
                   complex<T>(-0.6506315378609466, +0.7920349342629495),
                   complex<T>(-0.5564766488918566, +0.9148198405846728),
                   complex<T>(-0.4345168906815268, +1.045382255856986),
                   complex<T>(-0.2640041595834027, +1.192762031948052) },
                 1.0 };
    case 22:
        return { {},
                 { complex<T>(-0.2566376987939318, -1.1979824335552132),
                   complex<T>(-0.4232528745642629, -1.055755605227546),
                   complex<T>(-0.5430983056306306, -0.9299947824439877),
                   complex<T>(-0.6362427683267828, -0.811887504024635),
                   complex<T>(-0.7105305456418792, -0.6982266265924527),
                   complex<T>(-0.7700332930556816, -0.5874255426351151),
                   complex<T>(-0.8171682088462721, -0.47856194922027806),
                   complex<T>(-0.8534754036851689, -0.37103893194823206),
                   complex<T>(-0.8799661455640174, -0.26443630392015344),
                   complex<T>(-0.8972983138153532, -0.15843519122898653),
                   complex<T>(-0.9058702269930871, -0.05277490828999903),
                   complex<T>(-0.9058702269930871, +0.05277490828999903),
                   complex<T>(-0.8972983138153532, +0.15843519122898653),
                   complex<T>(-0.8799661455640174, +0.26443630392015344),
                   complex<T>(-0.8534754036851689, +0.37103893194823206),
                   complex<T>(-0.8171682088462721, +0.47856194922027806),
                   complex<T>(-0.7700332930556816, +0.5874255426351151),
                   complex<T>(-0.7105305456418792, +0.6982266265924527),
                   complex<T>(-0.6362427683267828, +0.811887504024635),
                   complex<T>(-0.5430983056306306, +0.9299947824439877),
                   complex<T>(-0.4232528745642629, +1.055755605227546),
                   complex<T>(-0.2566376987939318, +1.1979824335552132) },
                 1.0 };
    case 23:
        return { {},
                 { complex<T>(-0.24976972022089572, -1.2028131878706978),
                   complex<T>(-0.4126986617510148, -1.0653287944755134),
                   complex<T>(-0.5304922463810198, -0.9439760364018306),
                   complex<T>(-0.6225903228771341, -0.830155830281298),
                   complex<T>(-0.6965966033912708, -0.720734137475305),
                   complex<T>(-0.7564660146829886, -0.6141594859476034),
                   complex<T>(-0.8045561642053178, -0.5095305912227259),
                   complex<T>(-0.8423805948021129, -0.4062657948237603),
                   complex<T>(-0.8709469395587415, -0.3039581993950041),
                   complex<T>(-0.8909283242471254, -0.20230246993812237),
                   complex<T>(-0.9027564979912508, -0.10105343353140452),
                   complex<T>(-0.9066732476324991, +0.),
                   complex<T>(-0.9027564979912508, +0.10105343353140452),
                   complex<T>(-0.8909283242471254, +0.20230246993812237),
                   complex<T>(-0.8709469395587415, +0.3039581993950041),
                   complex<T>(-0.8423805948021129, +0.4062657948237603),
                   complex<T>(-0.8045561642053178, +0.5095305912227259),
                   complex<T>(-0.7564660146829886, +0.6141594859476034),
                   complex<T>(-0.6965966033912708, +0.720734137475305),
                   complex<T>(-0.6225903228771341, +0.830155830281298),
                   complex<T>(-0.5304922463810198, +0.9439760364018306),
                   complex<T>(-0.4126986617510148, +1.0653287944755134),
                   complex<T>(-0.24976972022089572, +1.2028131878706978) },
                 1.0 };
    case 24:
        return { {},
                 { complex<T>(-0.24334813375248746, -1.2072986837319728),
                   complex<T>(-0.4027853855197519, -1.0741951965186751),
                   complex<T>(-0.518591457482032, -0.9569048385259057),
                   complex<T>(-0.6096221567378332, -0.8470292433077199),
                   complex<T>(-0.6832565803536519, -0.7415032695091649),
                   complex<T>(-0.7433392285088533, -0.6388084216222569),
                   complex<T>(-0.7921695462343489, -0.5380628490968016),
                   complex<T>(-0.8312326466813242, -0.4386985933597306),
                   complex<T>(-0.8615278304016355, -0.34032021126186246),
                   complex<T>(-0.8837358034555707, -0.24263352344013836),
                   complex<T>(-0.8983105104397872, -0.14540561338736102),
                   complex<T>(-0.9055312363372773, -0.0484400665404787),
                   complex<T>(-0.9055312363372773, +0.0484400665404787),
                   complex<T>(-0.8983105104397872, +0.14540561338736102),
                   complex<T>(-0.8837358034555707, +0.24263352344013836),
                   complex<T>(-0.8615278304016355, +0.34032021126186246),
                   complex<T>(-0.8312326466813242, +0.4386985933597306),
                   complex<T>(-0.7921695462343489, +0.5380628490968016),
                   complex<T>(-0.7433392285088533, +0.6388084216222569),
                   complex<T>(-0.6832565803536519, +0.7415032695091649),
                   complex<T>(-0.6096221567378332, +0.8470292433077199),
                   complex<T>(-0.518591457482032, +0.9569048385259057),
                   complex<T>(-0.4027853855197519, +1.0741951965186751),
                   complex<T>(-0.24334813375248746, +1.2072986837319728) },
                 1.0 };
    default:
        return { {}, {}, 1.f };
    }
}

namespace internal
{

template <typename T>
KFR_FUNCTION zpk<T> bilinear(const zpk<T>& filter, identity<T> fs)
{
    const T fs2 = 2.0 * fs;
    zpk<T> result;
    result.z = (fs2 + filter.z) / (fs2 - filter.z);
    result.p = (fs2 + filter.p) / (fs2 - filter.p);
    result.z.resize(result.p.size(), complex<T>(-1));
    result.k = filter.k * kfr::real(product(fs2 - filter.z) / product(fs2 - filter.p));
    return result;
}

template <typename T>
struct zero_pole_pairs
{
    complex<T> p1, p2, z1, z2;
};

template <typename T>
KFR_FUNCTION vec<T, 3> zpk2tf_poly(const complex<T>& x, const complex<T>& y)
{
    return { T(1), -(x.real() + y.real()), x.real() * y.real() - x.imag() * y.imag() };
}

template <typename T>
KFR_FUNCTION biquad_params<T> zpk2tf(const zero_pole_pairs<T>& pairs, identity<T> k)
{
    vec<T, 3> zz = k * zpk2tf_poly(pairs.z1, pairs.z2);
    vec<T, 3> pp = zpk2tf_poly(pairs.p1, pairs.p2);
    // return { zz[0], zz[1], zz[2], pp[0], pp[1], pp[2] };
    return { pp[0], pp[1], pp[2], zz[0], zz[1], zz[2] };
}

template <typename T>
KFR_FUNCTION univector<complex<T>> cplxreal(const univector<complex<T>>& list)
{
    univector<complex<T>> x = list;
    std::sort(x.begin(), x.end(), [](const complex<T>& a, const complex<T>& b) { return a.real() < b.real(); });
    T tol                        = std::numeric_limits<T>::epsilon() * 100;
    univector<complex<T>> result = x;
    for (size_t i = result.size(); i > 1; i--)
    {
        if (!isreal(result[i - 1]) && !isreal(result[i - 2]))
        {
            if (abs(result[i - 1].real() - result[i - 2].real()) < tol &&
                abs(result[i - 1].imag() + result[i - 2].imag()) < tol)
            {
                result.erase(result.begin() + i - 1);
                result[i - 2].imag(abs(result[i - 2].imag()));
            }
        }
    }
    return result;
}

template <typename T>
KFR_FUNCTION size_t nearest_real_or_complex(const univector<complex<T>>& list, const complex<T>& val,
                                            bool mustbereal = true)
{
    univector<complex<T>> filtered;
    for (complex<T> v : list)
    {
        if (isreal(v) == mustbereal)
        {
            filtered.push_back(v);
        }
    }
    TESTO_ASSERT(!filtered.empty());
    if (filtered.empty())
        return std::numeric_limits<size_t>::max();

    size_t minidx = 0;
    T minval      = cabs(val - filtered[0]);
    for (size_t i = 1; i < list.size(); i++)
    {
        T newminval = cabs(val - filtered[i]);
        if (newminval < minval)
        {
            minval = newminval;
            minidx = i;
        }
    }
    return minidx;
}

template <typename T>
KFR_FUNCTION int countreal(const univector<complex<T>>& list)
{
    int nreal = 0;
    for (complex<T> c : list)
    {
        if (c.imag() == 0)
            nreal++;
    }
    return nreal;
}

template <typename T>
KFR_FUNCTION zpk<T> lp2lp_zpk(const zpk<T>& filter, identity<T> wo)
{
    zpk<T> result;
    result.z = wo * filter.z;
    result.p = wo * filter.p;
    result.k = filter.k * pow(wo, filter.p.size() - filter.z.size());
    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> lp2hp_zpk(const zpk<T>& filter, identity<T> wo)
{
    zpk<T> result;
    result.z = wo / filter.z;
    result.p = wo / filter.p;
    result.z.resize(result.p.size(), T(0));
    result.k = filter.k * kfr::real(product(-filter.z) / product(-filter.p));
    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> lp2bp_zpk(const zpk<T>& filter, identity<T> wo, identity<T> bw)
{
    zpk<T> lowpass;
    lowpass.z = bw * 0.5 * filter.z;
    lowpass.p = bw * 0.5 * filter.p;

    zpk<T> result;
    result.z = concatenate(lowpass.z + csqrt(csqr(lowpass.z) - wo * wo),
                           lowpass.z - csqrt(csqr(lowpass.z) - wo * wo));
    result.p = concatenate(lowpass.p + csqrt(csqr(lowpass.p) - wo * wo),
                           lowpass.p - csqrt(csqr(lowpass.p) - wo * wo));

    result.z.resize(result.z.size() + filter.p.size() - filter.z.size(), T(0));
    result.k = filter.k * pow(bw, filter.p.size() - filter.z.size());

    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> lp2bs_zpk(const zpk<T>& filter, identity<T> wo, identity<T> bw)
{
    zpk<T> highpass;
    highpass.z = (bw * 0.5) / filter.z;
    highpass.p = (bw * 0.5) / filter.p;

    zpk<T> result;
    result.z = concatenate(highpass.z + csqrt(csqr(highpass.z) - wo * wo),
                           highpass.z - csqrt(csqr(highpass.z) - wo * wo));
    result.p = concatenate(highpass.p + csqrt(csqr(highpass.p) - wo * wo),
                           highpass.p - csqrt(csqr(highpass.p) - wo * wo));

    result.z.resize(result.z.size() + filter.p.size() - filter.z.size(), complex<T>(0, +wo));
    result.z.resize(result.z.size() + filter.p.size() - filter.z.size(), complex<T>(0, -wo));
    result.k = filter.k * kfr::real(product(-filter.z) / product(-filter.p));

    return result;
}

template <typename T>
KFR_FUNCTION T warp_freq(T frequency, T fs)
{
    frequency = 2 * frequency / fs;
    fs        = 2.0;
    T warped  = 2 * fs * tan(c_pi<T> * frequency / fs);
    return warped;
}

} // namespace internal

template <typename T>
KFR_FUNCTION zpk<T> iir_lowpass(const zpk<T>& filter, identity<T> frequency, identity<T> fs = 2.0)
{
    T warped = internal::warp_freq(frequency, fs);

    zpk<T> result = filter;
    result        = internal::lp2lp_zpk(result, warped);
    result        = internal::bilinear(result, 2.0); // fs = 2.0
    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> iir_highpass(const zpk<T>& filter, identity<T> frequency, identity<T> fs = 2.0)
{
    T warped = internal::warp_freq(frequency, fs);

    zpk<T> result = filter;
    result        = internal::lp2hp_zpk(result, warped);
    result        = internal::bilinear(result, 2.0); // fs = 2.0
    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> iir_bandpass(const zpk<T>& filter, identity<T> lowfreq, identity<T> highfreq,
                                 identity<T> fs = 2.0)
{
    T warpedlow  = internal::warp_freq(lowfreq, fs);
    T warpedhigh = internal::warp_freq(highfreq, fs);

    zpk<T> result = filter;
    result        = internal::lp2bp_zpk(result, std::sqrt(warpedlow * warpedhigh), warpedhigh - warpedlow);
    result        = internal::bilinear(result, 2.0); // fs = 2.0
    return result;
}

template <typename T>
KFR_FUNCTION zpk<T> iir_bandstop(const zpk<T>& filter, identity<T> lowfreq, identity<T> highfreq,
                                 identity<T> fs = 2.0)
{
    T warpedlow  = internal::warp_freq(lowfreq, fs);
    T warpedhigh = internal::warp_freq(highfreq, fs);

    zpk<T> result = filter;
    result        = internal::lp2bs_zpk(result, std::sqrt(warpedlow * warpedhigh), warpedhigh - warpedlow);
    result        = internal::bilinear(result, 2.0); // fs = 2.0
    return result;
}

template <typename T>
KFR_FUNCTION std::vector<biquad_params<T>> to_sos(const zpk<T>& filter)
{
    if (filter.p.empty() && filter.z.empty())
        return { biquad_params<T>(filter.k, 0., 0., 1., 0., 0) };

    zpk<T> filt   = filter;
    size_t length = std::max(filter.p.size(), filter.z.size());
    filt.p.resize(length, complex<T>(0));
    filt.z.resize(length, complex<T>(0));

    size_t n_sections = (length + 1) / 2;
    if (length & 1)
    {
        filt.z.push_back(complex<T>(0));
        filt.p.push_back(complex<T>(0));
    }

    filt.z = internal::cplxreal(filt.z);
    filt.p = internal::cplxreal(filt.p);
    std::vector<internal::zero_pole_pairs<T>> pairs(n_sections);

    for (size_t si = 0; si < n_sections; si++)
    {
        size_t worstidx = 0;
        T worstval      = abs(1 - cabs(filt.p[0]));
        for (size_t i = 1; i < filt.p.size(); i++)
        {
            T val = abs(1 - cabs(filt.p[i]));
            if (val < worstval)
            {
                worstidx = i;
                worstval = val;
            }
        }
        complex<T> p1 = filt.p[worstidx];
        filt.p.erase(filt.p.begin() + worstidx);

        complex<T> z1, p2, z2;
        if (isreal(p1) && internal::countreal(filt.p) == 0)
        {
            size_t z1_idx = internal::nearest_real_or_complex(filt.z, p1, true);
            z1            = filt.z[z1_idx];
            filt.z.erase(filt.z.begin() + z1_idx);
            p2 = z2 = 0;
        }
        else
        {
            size_t z1_idx;
            if (!isreal(p1) && internal::countreal(filt.z) == 1)
            {
                z1_idx = internal::nearest_real_or_complex(filt.z, p1, false);
            }
            else
            {
                size_t minidx = 0;
                T minval      = cabs(p1 - filt.z[0]);
                for (size_t i = 1; i < filt.z.size(); i++)
                {
                    T newminval = cabs(p1 - filt.z[i]);
                    if (newminval < minval)
                    {
                        minidx = i;
                        minval = newminval;
                    }
                }
                z1_idx = minidx;
            }
            z1 = filt.z[z1_idx];
            filt.z.erase(filt.z.begin() + z1_idx);
            if (!isreal(p1))
            {
                if (!isreal(z1))
                {
                    p2 = cconj(p1);
                    z2 = cconj(z1);
                }
                else
                {
                    p2            = cconj(p1);
                    size_t z2_idx = internal::nearest_real_or_complex(filt.z, p1, true);
                    z2            = filt.z[z2_idx];
                    TESTO_ASSERT(isreal(z2));
                    filt.z.erase(filt.z.begin() + z2_idx);
                }
            }
            else
            {
                size_t p2_idx;
                size_t z2_idx;
                if (!isreal(z1))
                {
                    z2     = cconj(z1);
                    p2_idx = internal::nearest_real_or_complex(filt.z, p1, true);
                    p2     = filt.p[p2_idx];
                    TESTO_ASSERT(isreal(p2));
                }
                else
                {
                    size_t worstidx = 0;
                    T worstval      = abs(cabs(filt.p[0]) - 1);
                    for (size_t i = 1; i < filt.p.size(); i++)
                    {
                        T val = abs(cabs(filt.p[i]) - 1);
                        if (val < worstval)
                        {
                            worstidx = i;
                            worstval = val;
                        }
                    }
                    p2_idx = worstidx;
                    p2     = filt.p[p2_idx];

                    TESTO_ASSERT(isreal(p2));
                    z2_idx = internal::nearest_real_or_complex(filt.z, p2, true);
                    z2     = filt.z[z2_idx];
                    TESTO_ASSERT(isreal(z2));
                    filt.z.erase(filt.z.begin() + z2_idx);
                }
                filt.p.erase(filt.p.begin() + p2_idx);
            }
        }
        pairs[si].p1 = p1;
        pairs[si].p2 = p2;
        pairs[si].z1 = z1;
        pairs[si].z2 = z2;
    }

    std::vector<biquad_params<T>> result(n_sections);
    for (size_t si = 0; si < n_sections; si++)
    {
        result[si] = internal::zpk2tf(pairs[n_sections - 1 - si], si == 0 ? filt.k : T(1));
    }
    return result;
}
} // namespace CMT_ARCH_NAME

} // namespace kfr
