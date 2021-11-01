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
#include <cmath>		//std::isnan
#include <opencv2/core/core.hpp>
#include <Eigen/Core>

#include "AHCTypes.hpp"		//shared_ptr
#include "AHCParamSet.hpp"	//depthDisContinuous
#include "DisjointSet.hpp"	//PlaneSeg::mergeNbsFrom
#include "eigen33.h"					\

namespace ahc
{

//return true if d0 and d1 is discontinuous
template <class T> inline static bool
depthDisContinuous(T d0, T d1, const ParamSet<T>& params)
{
    return std::abs(d0 - d1) > params.T_dz(d0);
}

template <class T> inline static bool
is_valid(T z)
{
    return (z != 0 && !std::isnan(z));
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
    using vector3_t	= cv::Vec<T, 3>;
    using vector9_t	= cv::Vec<T, 9>;
  // using vector3_t	= Eigen::Matrix<T, 3, 1>;
  // using vector9_t	= Eigen::Matrix<T, 9, 1>;

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
		Stats()	:_sums(), _npoints(0)		{ clear(); }
		Stats(const Stats& a, const Stats& b)
		    :_sums(a._sums + b._sums),
		     _npoints(a._npoints + b._npoints)	{}


	size_t	npoints()			const	{ return _npoints; }

	void	clear()
		{
		    _sums.all(0);
		  // _sums.setZero();
		    _npoints = 0;
		}

	Stats&	push(const vector3_t& point)
		{
		    _sums += moment(point);
		    ++_npoints;

		    return *this;
		}
      /*
	Stats&	pop(const vector3_t& point)
		{
		    _sums -= moment(point);
		    --_npoints;

		    return *this;
		}

	Stats&	operator +=(const Stats& stats)
		{
		    _sums    += stats._sums;
		    _npoints += stats._npoints;

		    return *this;
		}

	Stats&	operator -=(const Stats& stats)
		{
		    _sums    -= stats._sums;
		    _npoints -= stats._npoints;

		    return *this;
		}
      */
	bool	compute(vector3_t& center,
			vector3_t& normal, T& mse, T& curvature) const
		{
		    assert(_npoints >= 4);

		    const T	sc = 1.0/_npoints;

		    center(0) = _sums(0)*sc;
		    center(1) = _sums(1)*sc;
		    center(2) = _sums(2)*sc;

		    T	K[3][3] = {{_sums(3) - _sums(0)*_sums(0)*sc,
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

		    normal(0) = V[0][0];
		    normal(1) = V[1][0];
		    normal(2) = V[2][0];

		  // enforce dot(normal,center)<00 so normal always
		  // points towards camera
		    if (normal.dot(center) > 0)
			normal *= -1;

		    mse       = evals[0] * sc;	// MSE (mean square error)
		    curvature = evals[0] / (evals[0] + evals[1] + evals[2]);

		    return true;
		}

      private:
	static vector9_t
		moment(const vector3_t& p)
		{
		    vector9_t	m;
		    m(0) = p(0);
		    m(1) = p(1);
		    m(2) = p(2);
		    m(3) = p(0) * p(0);
		    m(4) = p(0) * p(1);
		    m(5) = p(0) * p(2);
		    m(6) = p(1) * p(1);
		    m(7) = p(1) * p(2);
		    m(8) = p(2) * p(2);

		    return m;
		}

	static bool
		eig33sym(T K[3][3], T V[3][3], T evals[3])
		{
		    T	tmpV[3][3];
		    if (!eigen33(K, tmpV, evals) != 0)
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
	size_t		_npoints;
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
    PlaneSeg(const CLOUD& cloud, int root_block_id,
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
		const auto&	p = cloud.get(v, u);
		if (!is_valid(p(2)))
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
		if (u + 1 < imgWidth)
		{
		    const auto&	pn = cloud.get(v, u + 1);

		    if (is_valid(pn(2)) &&
			depthDisContinuous(p(2), pn(2), params))
		    {
			windowValid=false;
			break;
		    }
		}
		if (v + 1 < imgHeight)
		{
		    const auto& pn = cloud.get(v + 1, u);

		    if (is_valid(pn(2)) &&
			depthDisContinuous(p(2), pn(2), params))
		    {
			windowValid = false;
			break;
		    }
		}

		_stats.push(p);
	    }
	    if (!windowValid)
		break;
	}

	if (windowValid)
	{//if nan or depth-discontinuity shows, this obj will be rejected
	    _N	   = _stats.npoints();
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
	 _N(_stats.npoints()),
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
		    return std::abs(_normal.dot(p._normal));
		}

  /**
   *  \brief	signed distance between this plane and the point pt[3]
   */
    T		signedDist(const vector3_t pt) const
		{
		    return _normal.dot(pt - _center);
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
    vector3_t	_center; 	//q: plane center (center of mass)
    vector3_t	_normal; 	//n: plane equation n'p=q
    T		_mse;		//mean square error
    int		_N;		//#member points, same as _stats._N
    bool	_nouse;		//this PlaneSeg will be marked as _nouse after merged with others to produce a new PlaneSeg node in the graph

    std::set<Ptr>	_nbs;	//neighbors, i.e. adjacency list for a graph structure
};//PlaneSeg

}//ahc
