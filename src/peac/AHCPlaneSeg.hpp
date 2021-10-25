//
// Copyright 2014 Mitsubishi Electric Research Laboratories All
// Rights Reserved.
//
// Permission to use, copy and modify this software and its
// documentation without fee for educational, research and non-profit
// purposes, is hereby granted, provided that the above copyright
// notice, this paragraph, and the following three paragraphs appear
// in all copies.
//
// To request permission to incorporate this software into commercial
// products contact: Director; Mitsubishi Electric Research
// Laboratories (MERL); 201 Broadway; Cambridge, MA 02139.
//
// IN NO EVENT SHALL MERL BE LIABLE TO ANY PARTY FOR DIRECT,
// INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
// LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
// DOCUMENTATION, EVEN IF MERL HAS BEEN ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGES.
//
// MERL SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN
// "AS IS" BASIS, AND MERL HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
// SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
#pragma once

#include <set>			//PlaneSeg::NbSet
#include <vector>		//_mseseq
#include <limits>		//quiet_NaN

#include "AHCTypes.hpp"		//shared_ptr
#include "eig33sym.hpp"		//PlaneSeg::Stats::compute
#include "AHCParamSet.hpp"	//depthDisContinuous
#include "DisjointSet.hpp"	//PlaneSeg::mergeNbsFrom

namespace ahc
{

//return true if d0 and d1 is discontinuous
template <class T> inline static bool
depthDisContinuous(T d0, T d1, const ParamSet<T>& params)
{
    return std::abs(d0 - d1) > params.T_dz(d0);
}

/**
 *  \brief	PlaneSeg is a struct representing a Plane Segment
 *		as a node of a graph
 *  \details	It is usually dynamically allocated and garbage collected
 *		by boost::shared_ptr
 */
template <class T>
class PlaneSeg
{
  public:
    using Ptr		= PlaneSeg*;
    using shared_ptr	= ahc::shared_ptr<PlaneSeg>;
    using param_set_t	= ParamSet<T>;

  private:
  /**
   *  \brief	An internal struct holding this PlaneSeg's member points'
   *		1st and 2nd order statistics
   *  \details	It is usually dynamically allocated and garbage collected
   *		by boost::shared_ptr
   */
    class Stats
    {
      public:
	Stats()
	    :_sx(0), _sy(0), _sz(0),
	     _sxx(0), _syy(0), _szz(0), _sxy(0), _syz(0), _sxz(0),
	     _size(0)							{}

      //merge from two other Stats
	Stats(const Stats& a, const Stats& b)
	    :_sx(a._sx + b._sx), _sy(a._sy + b._sy), _sz(a._sz + b._sz),
	     _sxx(a._sxx + b._sxx), _syy(a._syy + b._syy),
	     _szz(a._szz + b._szz), _sxy(a._sxy + b._sxy),
	     _syz(a._syz + b._syz), _sxz(a._sxz + b._sxz),
	     _size(a._size + b._size)					{}

	int	size()				const	{ return _size; }
	void	clear()
		{
		    _sx = _sy = _sz
			= _sxx = _syy = _szz = _sxy = _syz = _sxz = 0;
		    _size = 0;
		}

      //push a new point (x,y,z) into this Stats
	auto&	push(T x, T y, T z)
		{
		    _sx += x;
		    _sy += y;
		    _sz += z;

		    _sxx += x*x;
		    _syy += y*y;
		    _szz += z*z;
		    _sxy += x*y;
		    _syz += y*z;
		    _sxz += x*z;

		    ++_size;

		    return *this;
		}

      //push a new Stats into this Stats
	auto&	operator +=(const Stats& other)
		{
		    _sx += other._sx;
		    _sy += other._sy;
		    _sz += other._sz;
		    _sxx += other._sxx;
		    _syy += other._syy;
		    _szz += other._szz;
		    _sxy += other._sxy;
		    _syz += other._syz;
		    _sxz += other._sxz;

		    _size += other._size;

		    return *this;
		}

      //caller is responsible to ensure (x,y,z) was collected in this stats
	auto&	pop(T x, T y, T z)
		{
		    _sx -= x;
		    _sy -= y;
		    _sz -= z;

		    _sxx -= x*x;
		    _syy -= y*y;
		    _szz -= z*z;
		    _sxy -= x*y;
		    _syz -= y*z;
		    _sxz -= x*z;

		    --_size;
		    assert(_size>=0);

		    return *this;
		}

      //caller is responsible to ensure {other} were collected in this stats
	auto&	operator -=(const Stats& other)
		{
		    _sx -= other._sx;
		    _sy -= other._sy;
		    _sz -= other._sz;

		    _sxx -= other._sxx;
		    _syy -= other._syy;
		    _szz -= other._szz;
		    _sxy -= other._sxy;
		    _syz -= other._syz;
		    _sxz -= other._sxz;

		    _size -= other._size;
		    assert(_size >= 0);

		    return *this;
		}

      /**
       *  \brief		PCA-based plane fitting
       *
       *  \param center		[out] center of mass of the PlaneSeg
       *  \param normal		[out] unit normal vector of the PlaneSeg
       *			(ensure normal.z>=0)
       *  \param mse		[out] mean-square-error of the plane fitting
       *  \param curvature	[out] defined as in pcl
       */
	void	compute(T center[3], T normal[3], T& mse, T& curvature) const
		{
		    assert(_size >= 4);

		    const T	sc = T(1.0)/_size;	//this->ids.size();
		  //calc plane equation: center, normal and mse
		    center[0] = _sx*sc;
		    center[1] = _sy*sc;
		    center[2] = _sz*sc;
		    T	K[3][3] =
			{{_sxx-_sx*_sx*sc, _sxy-_sx*_sy*sc, _sxz-_sx*_sz*sc},
			 {0,		   _syy-_sy*_sy*sc, _syz-_sy*_sz*sc},
			 {0,		   0,		    _szz-_sz*_sz*sc}};
		    K[1][0]=K[0][1]; K[2][0]=K[0][2]; K[2][1]=K[1][2];
		    T	sv[3]	= {0, 0, 0};
		    T	V[3][3] = {0};
		    LA::eig33sym(K, sv, V); //!!! first eval is the least one
		  //LA.svd33(K, sv, V);
		    if (V[0][0]*center[0] +
			V[1][0]*center[1] + V[2][0]*center[2] <= 0)
		    {
		      // enforce dot(normal,center)<00 so normal always
		      // points towards camera
			normal[0] = V[0][0];
			normal[1] = V[1][0];
			normal[2] = V[2][0];
		    }
		    else
		    {
			normal[0] = -V[0][0];
			normal[1] = -V[1][0];
			normal[2] = -V[2][0];
		    }
		    mse	      = sv[0] * sc;
		    curvature = sv[0]/(sv[0] + sv[1] + sv[2]);
		}

      private:
	T	_sx, _sy, _sz,		// sum of x/y/z
		_sxx, _syy, _szz,	// sum of xx/yy/zz
		_sxy, _syz, _sxz;	// sum of xy/yz/xz
	int	_size;			// #points in this PlaneSeg
    };

  public:
  /**
   *  \brief construct a PlaneSeg during graph initialization
   *
   *  \param points		[in] organized point cloud adapter,
   *				see NullImage3D
   *  \param root_block_id	[in] initial window/block's id
   *  \param seed_row		[in] row index of the upper left pixel
   *				of the initial window/block
   *  \param seed_col		[in] row index of the upper left pixel
   *				of the initial window/block
   *  \param imgWidth		[in] width of the organized point cloud
   *  \param imgHeight		[in] height of the organized point cloud
   *  \param winWidth		[in] width of the initial window/block
   *  \param winHeight		[in] height of the initial window/block
   *  \param depthChangeFactor	[in] parameter to determine depth discontinuity
   *
   *  \details if exist depth discontinuity in this initial PlaneSeg, nouse will be set true and N 0.
   */
    template<class CLOUD>
    PlaneSeg(const CLOUD& points, int root_block_id,
	     int seed_row, int seed_col, int imgWidth, int imgHeight,
	     int winWidth, int winHeight, const param_set_t& params)
	:_stats(),
	 _curvature(0),
	 _rid(root_block_id),
	 _center(),
	 _normal(),
	 _mse(),
	 _N(0),
	 _nouse(false),
	 _nbs()
    {
      //assert(0<=seed_row && seed_row<height && 0<=seed_col && seed_col<width && winW>0 && winH>0);
	bool		windowValid = true;
	int		nanCnt = 0, nanCntTh = winHeight*winWidth/2;
      //calc _stats
	const auto	vmax = std::min(imgHeight, seed_row + winHeight);
	const auto	umax = std::min(imgWidth,  seed_col + winWidth);

	for (int v = seed_row; v < vmax; ++v)
	{
	    for (int u = seed_col; u < umax; ++u)
	    {
		T	x = 0, y = 0, z = 10000;
		if (!points.get(v, u, x, y, z))
		{
		    if (params.initType() == param_set_t::INIT_LOOSE)
		    {
			++nanCnt;
			if (nanCnt < nanCntTh)
			    continue;
		    }
		    windowValid=false;
		    break;
		}
		T	xn = 0, yn = 0, zn = 10000;
		if (u + 1 < imgWidth &&
		    points.get(v, u + 1, xn, yn, zn) &&
		    depthDisContinuous(z, zn, params))
		{
		    windowValid=false;
		    break;
		}
		if (v + 1 < imgHeight &&
		    points.get(v + 1, u, xn, yn, zn) &&
		    depthDisContinuous(z, zn, params))
		{
		    windowValid = false;
		    break;
		}
		_stats.push(x, y, z);
	    }
	    if (!windowValid)
		break;
	}

	if (windowValid)
	{//if nan or depth-discontinuity shows, this obj will be rejected
	    _N	   = _stats.size();
	    _nouse = false;
	}
	else
	{
	    _stats.clear();
	    _N = 0;
	    _nouse = true;
	}

	if (_N < 4)
	{
	    _mse = _curvature = std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
	    _stats.compute(_center, _normal, _mse, _curvature);
	  //nbs information to be maintained outside the class
	  //typically when initializing the graph structure
	}
      //std::cout<<_curvature<<std::endl;
    }

  /**
   *  \brief	construct a new PlaneSeg from two PlaneSeg pa and pb
   *		when trying to merge
   *
   *  \param pa	[in] a PlaneSeg
   *  \param pb	[in] a PlaneSeg
   */
    PlaneSeg(const PlaneSeg& pa, const PlaneSeg& pb)
	:_stats(pa._stats, pb._stats),
	 _rid(pa._N >= pb._N ? pa._rid : pb._rid),
	 _N(_stats.size()),
	 _nouse(false)
    {
      //ds.union(pa._rid, pb._rid) will be called later
      //in mergeNbsFrom(pa,pb) function, since
      //this object might not be accepted into the graph structure

	_stats.compute(_center, _normal, _mse, _curvature);

      //nbs information to be maintained later if this node is accepted
    }

    void	update()
		{
		    _stats.compute(_center, _normal,_mse, _curvature);
		}

  /**
   *  \brief	similarity of two plane normals
   *
   *  \param p	[in] another PlaneSeg
   *  \return	abs(dot(_normal, p->_normal))
   *
   *  \details	1 means identical, 0 means perpendicular
   */
    T		normalSimilarity(const PlaneSeg& p) const
		{
		    return std::abs(_normal[0] * p._normal[0] +
				    _normal[1] * p._normal[1] +
				    _normal[2] * p._normal[2]);
		}

  /**
   *  \brief	signed distance between this plane and the point pt[3]
   */
    T		signedDist(const T pt[3]) const
		{
		    return _normal[0] * (pt[0] - _center[0]) +
			   _normal[1] * (pt[1] - _center[1]) +
			   _normal[2] * (pt[2] - _center[2]);
		}

  /**
   *  \brief	connect this PlaneSeg to another PlaneSeg p in the graph
   *
   *  \param p	[in] the other PlaneSeg
   */
    void	connect(Ptr p)
		{
		    if (p)
		    {
			_nbs.insert(p);
			p->_nbs.insert(this);
		    }
		}

  /**
   *  \brief	disconnect this PlaneSeg with all its neighbors
   *
   *  \details	after this call, _nbs.nbs should not contain this,
   *		and _nbs should be empty i.e. after this call
   *		this PlaneSeg node should be isolated in the graph
   */
    void	disconnectAllNbs()
		{
		    for (auto nb : _nbs)
			if (!nb->_nbs.erase(this))
			    std::cout << "[PlaneSeg warn] _nbs._nbs"
				" should have contained this!"
				      << std::endl;
		    _nbs.clear();
		}

  /**
   *  \brief finish merging PlaneSeg pa and pb to this
   *
   *  \param pa	[in] a parent PlaneSeg of this
   *  \param pb	[in] another parent PlaneSeg of this
   *  \param ds	[in] the disjoint set of initial window/block membership
   *		to be updated
   *
   *  \details	Only call this if this obj is accepted to be added
   *		to the graph of PlaneSeg pa and pb should not exist
   *		after this function is called, i.e. after this call
   *		this PlaneSeg node will be representing a merged node
   *		of pa and pb, and pa/pb will be isolated (and thus
   *		Garbage Collected) in the graph
   */
    void	mergeNbsFrom(PlaneSeg& pa, PlaneSeg& pb, DisjointSet& ds)
		{
		  //now we are sure that merging pa and pb is accepted
		    ds.Union(pa._rid, pb._rid);

		  //the new neighbors should be pa._nbs+pb._nbs-pa-pb
		    _nbs.insert(pa._nbs.begin(), pa._nbs.end());
		    _nbs.insert(pb._nbs.begin(), pb._nbs.end());
		    _nbs.erase(&pa);
		    _nbs.erase(&pb);

		  //pa and pb should be GC later after the following two steps
		    pa.disconnectAllNbs();
		    pb.disconnectAllNbs();

		  //complete the neighborhood from the other side
		    for (auto nb : _nbs)
			nb->_nbs.insert(this);

		    pa._nouse=pb._nouse=true;
		}

  private:
    Stats	_stats;
    T		_curvature;

  public:
    int		_rid;		//root block id
    T		_center[3]; 	//q: plane center (center of mass)
    T		_normal[3]; 	//n: plane equation n'p=q
    T		_mse;		//mean square error
    int		_N;		//#member points, same as _stats._N
    bool	_nouse;		//this PlaneSeg will be marked as _nouse after merged with others to produce a new PlaneSeg node in the graph

    std::set<Ptr>	_nbs;	//neighbors, i.e. adjacency list for a graph structure
};//PlaneSeg

}//ahc
