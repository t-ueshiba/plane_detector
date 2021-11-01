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
tridiagonal33(const T A[3][3], T Q[3][3], T d[3], T e[2])
{
  // -----------------------------------------------------------------------
  // Reduces a symmetric 3x3 matrix to tridiagonal form by applying
  // (unitary) Householder transformations:
  //            [ d[0]  e[0]       ]
  //    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
  //            [       e[1]  d[2] ]
  // The function accesses only the diagonal and upper triangular parts of
  // A. The access is read-only.
  // -----------------------------------------------------------------------

  // Initialize Q to the identitity matrix
    for (int i = 0; i < 3; ++i)
    {
	Q[i][i] = 1.0;
	for (int j=0; j < i; ++j)
	    Q[i][j] = Q[j][i] = 0.0;
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

      // Store inverse Householder transformation in Q
	for (int j = 1; j < 3; ++j)
	{
	    f = omega * u[j];
	    for (int i = 1; i < 3; ++i)
		Q[i][j] = Q[i][j] - f*u[i];
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
qr33(const T A[3][3], T Q[3][3], T w[3])
{
  // Transform A to real tridiagonal form by the Householder method
    T e[3];                   // The third element is used only as temporary
    tridiagonal33(A, Q, w, e);

  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
    for (int l = 0; l < 2; ++l)
    {
	for (int nIter = 0; ; )
	{
	    int	m;
	  // Check for convergence and exit iteration loop if off-diagonal
	  // element e(l) is zero
	    for (m = l; m <= 1; ++m)
	    {
		const T	g = std::abs(w[m])+std::abs(w[m+1]);
		if (std::abs(e[m]) + g == g)
		    break;
	    }
	    if (m == l)
		break;

	    if (nIter++ >= 30)
		return false;

	  // Calculate g = d_m - k
	    T	g = (w[l+1] - w[l]) / (e[l] + e[l]);
	    T	r = std::sqrt(square(g) + 1.0);
	    if (g > 0)
		g = w[m] - w[l] + e[l]/(g + r);
	    else
		g = w[m] - w[l] + e[l]/(g - r);

	    T	s = 1.0;
	    T	c = 1.0;
	    T	p = 0.0;
	    for (int i = m - 1; i >= l; --i)
	    {
		const T	f = s * e[i];
		const T	b = c * e[i];
		if (std::abs(f) > std::abs(g))
		{
		    c      = g / f;
		    r      = std::sqrt(square(c) + 1.0);
		    e[i+1] = f * r;
		    c     *= (s = 1.0/r);
		}
		else
		{
		    s      = f / g;
		    r      = std::sqrt(square(s) + 1.0);
		    e[i+1] = g * r;
		    s     *= (c = 1.0/r);
		}

		g = w[i+1] - p;
		r = (w[i] - g)*s + 2.0*c*b;
		p = s * r;
		w[i+1] = g + p;
		g = c*r - b;

	      // Form eigenvectors
		for (int k = 0; k < 3; ++k)
		{
		    const auto	t = Q[k][i+1];
		    Q[k][i+1] = s*Q[k][i] + c*t;
		    Q[k][i]   = c*Q[k][i] - s*t;
		}
	    }
	    w[l] -= p;
	    e[l]  = g;
	    e[m]  = 0.0;
	}
    }

    return true;
}

template <class T> bool
eigen33(const T A[3][3], T Q[3][3], T w[3])
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

    Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
    Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0];
    Q[2][1] = square(A[0][1]);

  // Calculate first eigenvector by the formula
  //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
    Q[0][0] = Q[0][1] + A[0][2]*w[0];
    Q[1][0] = Q[1][1] + A[1][2]*w[0];
    Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
    T	norm = square(Q[0][0]) + square(Q[1][0]) + square(Q[2][0]);

  // If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A[i][i] - w[0], fall
  // back to QL algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
  // i.e. all columns of A - w[0] * I are linearly dependent.
    if (norm <= error)
	return qr33(A, Q, w);
    else                      // This is the standard branch
    {
	norm = std::sqrt(1.0 / norm);
	for (int j=0; j < 3; ++j)
	    Q[j][0] *= norm;
    }

  // Calculate second eigenvector by the formula
  //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
    Q[0][1] = Q[0][1] + A[0][2]*w[1];
    Q[1][1] = Q[1][1] + A[1][2]*w[1];
    Q[2][1] = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
    norm    = square(Q[0][1]) + square(Q[1][1]) + square(Q[2][1]);
    if (norm <= error)
	return qr33(A, Q, w);
    else
    {
	norm = std::sqrt(1.0 / norm);
	for (int j = 0; j < 3; ++j)
	    Q[j][1] *= norm;
    }

  // Calculate third eigenvector according to
  //   v[2] = v[0] x v[1]
    Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1];
    Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1];
    Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1];

    return true;
}
