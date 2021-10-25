/*!
  \file		Region.h
  \author	Toshio UESHIBA
*/
#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace plane_detector
{
/************************************************************************
*  class Region								*
************************************************************************/
template <class T>
class Region
{
  public:
    using vector3_t	= Eigen::Matrix<T, 3, 1>;
    using vector9_t	= Eigen::Matrix<T, 9, 1>;

  public:
    template <class ITER>
		Region(ITER ps, ITER pe)
		    :_sums(vector9_t::Zero()), _npoints(0),
		     _center(vector3_t::Zero()),
		     _normal(vector3_t::Zero()), _mse(0)
		{
		    for (auto p = ps; p != pe; ++p)
			add(*p);
		    fit_plane();
		}
    
    size_t	npoints()			const	{ return _npoints; }
    T		mse()				const	{ return _mse; }

    Region	merge(const Region& region)
		{
		    auto	new_region = *this;
		    new_region += region.region;
		    new_region.fit_plane();

		    return new_region;
		}

  private:
    void	clear()
		{
		    _sums.setZero();
		    _npoints = 0;
		    _center.setZero();
		    _normal.setZero();
		    _mse = 0;
		}

    Region&	add(const vector3_t& point)
		{
		    _sums += moment(point);
		    ++_npoints;

		    return *this;
		}
    
    Region&	sub(const vector3_t& point)
		{
		    _sums -= moment(point);
		    --_npoints;

		    return *this;
		}

    Region&	operator +=(const Region& region)
		{
		    _sums    += region._sums;
		    _npoints += region._npoints;

		    return *this;
		}
    
    Region&	operator -=(const Region& region)
		{
		    _sums    -= region._sums;
		    _npoints -= region._npoints;

		    return *this;
		}

    bool	fit_plane() const
		{
		    assert(_size >= 4);

		    const T	sc = 1.0/_size;

		    _center = _sums.segment<3>(0)* sc;

		    const T	K[3][3] = {{_sums(3) - _sums(0)*_sums(0)*sc,
					    _sums(4) - _sums(0)*_sums(1)*sc,
					    _sums(5) - _sums(0)*_sums(2)*sc},
					   {0,
					    _sums(6) - _sums(1)*_sums(1)*sc,
					    _sums(7) - _sums(1)*_sums(2)*sc},
					   {0,
					    0,
					    _sums(8) - _sums(2)*_sums(2)*sc}};
		    // K[1][0] = K[0][1];
		    // K[2][0] = K[0][2];
		    // K[2][1] = K[1][2];
		    T	V[3][3], evals[3];
		    if (!eig33sym(K, V, evals))	// first eval is the least one
			return false;

		    _normal(0) = V[0][0];
		    _normal(1) = V[1][0];
		    _normal(2) = V[2][0];
		    
		  // enforce dot(normal,center)<00 so normal always
		  // points towards camera
		    if (_normal.dot(_center) > 0)
			_normal *= -1;

		    _mse       = evals[0] * sc;	// MSE (mean square error)
		    _curvature = evals[0] / (evals[0] + evals[1] + evals[2]);

		    return true;
		}

    static vector9_t
		moment(const vector3_t& p)
		{
		    vector9_t	m;
		    m.segment<3>(0) = p;
		    m.segment<3>(3) = p(0) * p;
		    m.segment<2>(6) = p(1) * p.segment<2>(1);
		    m(8)	    = p(2) * p(2);
		    
		    return m;
		}

    static bool	eig33sym(const T K[3][3], T V[3][3], T evals[3])
		{
		    T	tmpV[3][3];
		    if (dsyevh3(K, tmpV, evals) != 0)
			return false;

		    int	order[] = {0, 1, 2};
		    for(int i = 0; i < 3; ++i)
			for(int j = i + 1; j < 3; ++j)
			    if (evals[i] > evals[j])
			    {
				std::swap(evals[i], evals[j]);
				std::swap(order[i], order[j]);
			    }
		    V[0][0] = tmpV[0][order[0]];
		    V[0][1] = tmpV[0][order[1]];
		    V[0][2] = tmpV[0][order[2]];
		    V[1][0] = tmpV[1][order[0]];
		    V[1][1] = tmpV[1][order[1]];
		    V[1][2] = tmpV[1][order[2]];
		    V[2][0] = tmpV[2][order[0]];
		    V[2][1] = tmpV[2][order[1]];
		    V[2][2] = tmpV[2][order[2]];

		    return true;
		}
    
  private:
    vector9_t	_sums;
    size_t	_npoints;
    vector3_t	_center;	// q: plane center (center of mass)
    vector3_t	_normal;	// n: plane equation n'p=q
    T		_mse;		// mean square error
    T		_curvature;
};

}