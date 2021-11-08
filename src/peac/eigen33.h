// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#pragma once

#include <cmath>
#include <limits>

#define M_SQRT3	1.73205080756887729352744634151		// sqrt(3)

template <class T> static inline T	square(T x) { return x*x; }

template <class T> void
cardano(const T A[3][3], T w[3])
{
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
    const T	de = A[0][1] * A[1][2];		// d * e
    const T	dd = square(A[0][1]);		// d^2
    const T	ee = square(A[1][2]);		// e^2
    const T	ff = square(A[0][2]);		// f^2
    const T	m  = A[0][0] + A[1][1] + A[2][2];
    const T	c1 = (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2])
		   - (dd + ee + ff);
    const T	c0 = A[2][2]*dd + A[0][0]*ee + A[1][1]*ff
		   - A[0][0]*A[1][1]*A[2][2] - 2.0 * A[0][2]*de;

    const T	p      = square(m) - 3.0*c1;
    const T	q      = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
    const T	sqrt_p = std::sqrt(std::abs(p));

    T	phi = 27.0 * (0.25*square(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
    phi = (1.0/3.0) * std::atan2(std::sqrt(std::abs(phi)), q);

    const T	c = sqrt_p*std::cos(phi);
    const T	s = (1.0/M_SQRT3)*sqrt_p*std::sin(phi);

    w[1]  = (1.0/3.0)*(m - c);
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;
}

template <class T> void
tridiagonal33(const T A[3][3], T Qt[3][3], T d[3], T e[2])
{
  // -----------------------------------------------------------------------
  // Reduces a symmetric 3x3 matrix to tridiagonal form by applying
  // (unitary) Householder transformations:
  //             [ d[0]  e[0]       ]
  //    A = Qt . [ e[0]  d[1]  e[1] ] . Q
  //             [       e[1]  d[2] ]
  // The function accesses only the diagonal and upper triangular parts of
  // A. The access is read-only.
  // -----------------------------------------------------------------------

  // Initialize Qt to the identitity matrix
    for (int i = 0; i < 3; ++i)
    {
	Qt[i][i] = 1.0;
	for (int j = 0; j < i; ++j)
	    Qt[i][j] = Qt[j][i] = 0.0;
    }

  // Bring first row and column to the desired form
    const T	h = square(A[0][1]) + square(A[0][2]);
    const T	g = (A[0][1] > 0 ? -std::sqrt(h) : std::sqrt(h));
    e[0] = g;

    T		f = g * A[0][1];
    T		omega = h - f;
    
    if (omega > 0.0)
    {
	T	u[3], q[3];
	u[1] = A[0][1] - g;
	u[2] = A[0][2];

	omega = 1.0 / omega;
	T	K = 0.0;
	for (int i = 1; i < 3; ++i)
	{
	    f    = A[1][i] * u[1] + A[i][2] * u[2];
	    q[i] = omega * f;                  // p
	    K   += u[i] * f;                   // u* A u
	}
	K *= 0.5 * square(omega);

	for (int i = 1; i < 3; ++i)
	    q[i] = q[i] - K * u[i];

	d[0] = A[0][0];
	d[1] = A[1][1] - 2.0*q[1]*u[1];
	d[2] = A[2][2] - 2.0*q[2]*u[2];

      // Store inverse Householder transformation in Qt
	for (int i = 1; i < 3; ++i)
	{
	    f = omega * u[i];
	    for (int j = 1; j < 3; ++j)
		Qt[i][j] = Qt[i][j] - f*u[j];
	}

      // Calculate updated A[1][2] and store it in e[1]
	e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
    }
    else
    {
	for (int i = 0; i < 3; ++i)
	    d[i] = A[i][i];
	e[1] = A[1][2];
    }
}

template <class T> bool
qr33(const T A[3][3], T Qt[3][3], T w[3])
{
  // Transform A to real tridiagonal form by the Householder method
    T e[3];                   // The third element is used only as temporary
    tridiagonal33(A, Qt, w, e);

  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
    for (int n = 0; n < 2; ++n)	// l = 0, 1
    {
	for (int nIter = 0; ; )
	{
	    int	i = n;
	  // Check for convergence and exit iteration loop if off-diagonal
	  // element e(n) is zero
	    for (; i < 2; ++i)
	    {
		const auto	g = std::abs(w[i]) + std::abs(w[i+1]);
		if (std::abs(e[i]) + g == g)
		{
		    e[i] = 0;
		    break;
		}
	    }
	    if (i == n)
		break;

	    if (nIter++ >= 30)
		return false;

	  // [n = 0]: i = 1, 2; [n = 1]: i = 2	    
	    const auto	t = (w[n+1] - w[n])/(e[n] + e[n]);
	    const auto	r = std::sqrt(square(t) + T(1));
	    e[i] = w[i] - w[n] + e[n]/(t + (t > 0 ? r : -r));

	    auto	s = T(1);
	    auto	c = T(1);
	  // [n = 0, i = 1]: i = 0; [n = 0, i = 2]: i = 1, 0
	  // [n = 1, i = 2]: i = 1
	    while (--i >= n)
	    {
		const auto	x = e[i+1];
		const auto	y = s*e[i];
		const auto	z = c*e[i];
		if (std::abs(x) > std::abs(y))
		{
		    const auto	t = y/x;
		    const auto	r = std::sqrt(square(t) + T(1));
		    e[i+1] = x*r;
		    c	   = T(1)/r;
		    s      = c*t;
		}
		else
		{
		    const auto	t = x/y;
		    const auto	r = std::sqrt(square(t) + T(1));
		    e[i+1] = y*r;
		    s	   = T(1)/r;
		    c	   = s*t;
		}

		const auto	v = s*(w[i] - w[i+1]) + T(2)*c*z;
		const auto	p = s*v;
		w[i]   -= p;
		w[i+1] += p;
		e[i]    = c*v - z;	// updated e[i]
		    
	      // Form eigenvectors
		for (int k = 0; k < 3; ++k)
		{
		    const auto	t = Qt[i][k];
		    Qt[i][k]   = c*t - s*Qt[i+1][k];
		    Qt[i+1][k] = s*t + c*Qt[i+1][k];
		}
	    }
	}
    }

    return true;
}

template <class T> bool
eigen33(const T A[3][3], T Qt[3][3], T w[3])
{
  // Calculate eigenvalues
    cardano(A, w);

  //  n0 = square(A[0][0]) + square(A[0][1]) + square(A[0][2]);
  //  n1 = square(A[0][1]) + square(A[1][1]) + square(A[1][2]);

    T	t = std::abs(w[0]), u;
    if ((u=std::abs(w[1])) > t)
	t = u;
    if ((u=std::abs(w[2])) > t)
	t = u;
    if (t < 1.0)
	u = t;
    else
	u = square(t);
    const T	error = 256.0 * std::numeric_limits<T>::epsilon() * square(u);
  //  error = 256.0 * DBL_EPSILON * (n0 + u) * (n1 + u);

    Qt[1][0] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    Qt[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
    Qt[1][2] = square(A[0][1]);

  // Calculate first eigenvector by the formula
  //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
    Qt[0][0] = Qt[1][0] + A[0][2]*w[0];
    Qt[0][1] = Qt[1][1] + A[1][2]*w[0];
    Qt[0][2] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Qt[1][2];
    T	norm = square(Qt[0][0]) + square(Qt[0][1]) + square(Qt[0][2]);

  // If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A[i][i] - w[0], fall
  // back to QL algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
  // i.e. all columns of A - w[0] * I are linearly dependent.
    if (norm <= error)
    	return qr33(A, Qt, w);
    else                      // This is the standard branch
    {
	norm = std::sqrt(1.0 / norm);
	for (int j = 0; j < 3; ++j)
	    Qt[0][j] *= norm;
    }

  // Calculate second eigenvector by the formula
  //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
    Qt[1][0] = Qt[1][0] + A[0][2]*w[1];
    Qt[1][1] = Qt[1][1] + A[1][2]*w[1];
    Qt[1][2] = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Qt[1][2];
    norm    = square(Qt[1][0]) + square(Qt[1][1]) + square(Qt[1][2]);
    if (norm <= error)
  	return qr33(A, Qt, w);
    else
    {
	norm = std::sqrt(1.0 / norm);
	for (int j = 0; j < 3; ++j)
	    Qt[1][j] *= norm;
    }

  // Calculate third eigenvector according to
  //   v[2] = v[0] x v[1]
    Qt[2][0] = Qt[0][1]*Qt[1][2] - Qt[0][2]*Qt[1][1];
    Qt[2][1] = Qt[0][2]*Qt[1][0] - Qt[0][0]*Qt[1][2];
    Qt[2][2] = Qt[0][0]*Qt[1][1] - Qt[0][1]*Qt[1][0];

    return true;
}
