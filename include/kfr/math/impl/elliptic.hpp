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

namespace kfr
{
inline namespace CMT_ARCH_NAME
{

template <typename T>
KFR_INTRINSIC complex<T> sqrt(complex<T> w);
template <typename T>
KFR_INTRINSIC complex<T> atanh(complex<T> w);
template <typename T>
KFR_INTRINSIC complex<T> asin(complex<T> w);
template <typename T>
KFR_INTRINSIC T cmpl(T kx);
template <typename T>
KFR_INTRINSIC complex<T> cmpl(complex<T> kx);
template <typename T>
KFR_INTRINSIC T cephes_polevl(T x, T *coef, int N);

template <typename T>
KFR_INTRINSIC void sncndn ( T u, T m, T &sn, T &cn, T &dn );

template <typename T>
KFR_INTRINSIC complex<T> inv_jacobi_sn(complex<T> w, T m);
template <typename T>
KFR_INTRINSIC T inv_jacobi_sc1(T w, T m);

template <typename T>
KFR_INTRINSIC T elliptic_k (T m);
template <typename T>
KFR_INTRINSIC T elliptic_km1 (T p);
template <typename T>
KFR_INTRINSIC T elliptic_deg ( int N, T m1);

} // namespace CMT_ARCH_NAME
} // namespace kfr
