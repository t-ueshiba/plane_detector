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
    struct Stats
    {
	T	sx, sy, sz,	// sum of x/y/z
		sxx, syy, szz,	// sum of xx/yy/zz
		sxy, syz, sxz;	// sum of xy/yz/xz
	int	N;		// #points in this PlaneSeg

	Stats()
	    :sx(0), sy(0), sz(0),
	     sxx(0), syy(0), szz(0), sxy(0), syz(0), sxz(0), N(0)	{}

      //merge from two other Stats
	Stats(const Stats& a, const Stats& b)
	    :sx(a.sx+b.sx), sy(a.sy+b.sy), sz(a.sz+b.sz),
	     sxx(a.sxx+b.sxx), syy(a.syy+b.syy), szz(a.szz+b.szz),
	     sxy(a.sxy+b.sxy), syz(a.syz+b.syz), sxz(a.sxz+b.sxz),
	     N(a.N+b.N)							{}

	void	clear()
		{
		    sx = sy = sz = sxx = syy = szz = sxy = syz = sxz = 0;
		    N = 0;
		}

      //push a new point (x,y,z) into this Stats
	auto&	push(T x, T y, T z)
		{
		    sx += x;
		    sy += y;
		    sz += z;

		    sxx += x*x;
		    syy += y*y;
		    szz += z*z;
		    sxy += x*y;
		    syz += y*z;
		    sxz += x*z;

		    ++N;

		    return *this;
		}

      //push a new Stats into this Stats
	auto&	operator +=(const Stats& other)
		{
		    sx += other.sx;
		    sy += other.sy;
		    sz += other.sz;
		    sxx += other.sxx;
		    syy += other.syy;
		    szz += other.szz;
		    sxy += other.sxy;
		    syz += other.syz;
		    sxz += other.sxz;

		    N += other.N;

		    return *this;
		}

      //caller is responsible to ensure (x,y,z) was collected in this stats
	auto&	pop(T x, T y, T z)
		{
		    sx -= x;
		    sy -= y;
		    sz -= z;

		    sxx -= x*x;
		    syy -= y*y;
		    szz -= z*z;
		    sxy -= x*y;
		    syz -= y*z;
		    sxz -= x*z;

		    --N;
		    assert(N>=0);

		    return *this;
		}

      //caller is responsible to ensure {other} were collected in this stats
	auto&	operator -=(const Stats& other)
		{
		    sx -= other.sx;
		    sy -= other.sy;
		    sz -= other.sz;

		    sxx -= other.sxx;
		    syy -= other.syy;
		    szz -= other.szz;
		    sxy -= other.sxy;
		    syz -= other.syz;
		    sxz -= other.sxz;

		    N -= other.N;
		    assert(N>=0);

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
		    assert(N>=4);

		    const T	sc=T(1.0)/this->N;	//this->ids.size();
	  //calc plane equation: center, normal and mse
		    center[0] = sx*sc;
		    center[1] = sy*sc;
		    center[2] = sz*sc;
		    T	K[3][3] = {{sxx-sx*sx*sc, sxy-sx*sy*sc, sxz-sx*sz*sc},
				   {0,		  syy-sy*sy*sc, syz-sy*sz*sc},
				   {0,		  0,		szz-sz*sz*sc}};
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
		    mse = sv[0] * sc;
		    curvature = sv[0]/(sv[0] + sv[1] + sv[2]);
		}
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
    {
      //assert(0<=seed_row && seed_row<height && 0<=seed_col && seed_col<width && winW>0 && winH>0);
	this->rid = root_block_id;

	bool	windowValid = true;
	int	nanCnt = 0, nanCntTh = winHeight*winWidth/2;
      //calc _stats
	for (int i=seed_row, icnt=0; icnt < winHeight && i < imgHeight;
	     ++i, ++icnt)
	{
	    for (int j=seed_col, jcnt=0; jcnt < winWidth && j < imgWidth;
		 ++j, ++jcnt)
	    {
		T	x = 0, y = 0, z = 10000;
		if (!points.get(i, j, x, y, z))
		{
		    if (params.initType() == param_set_t::INIT_LOOSE)
		    {
			++nanCnt;
			if (nanCnt < nanCntTh)
			    continue;
		    }
#ifdef DEBUG_INIT
		    _type = TYPE_MISSING_DATA;
#endif
		    windowValid=false;
		    break;
		}
		T	xn = 0, yn = 0, zn = 10000;
		if (j + 1 < imgWidth &&
		    points.get(i, j+1, xn, yn, zn) &&
		    depthDisContinuous(z, zn, params))
		{
#ifdef DEBUG_INIT
		    _type = TYPE_DEPTH_DISCONTINUE;
#endif
		    windowValid=false;
		    break;
		}
		if (i + 1 < imgHeight &&
		    points.get(i + 1, j, xn, yn, zn) &&
		    depthDisContinuous(z,zn,params))
		{
#ifdef DEBUG_INIT
		    _type = TYPE_DEPTH_DISCONTINUE;
#endif
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
	    this->nouse=false;
	    this->N=_stats.N;
#ifdef DEBUG_INIT
	    _type = TYPE_NORMAL;
#endif
	}
	else
	{
	    this->N=0;
	    _stats.clear();
	    this->nouse=true;
	}

	if (this->N<4)
	{
	    this->mse = _curvature = std::numeric_limits<T>::quiet_NaN();
	}
	else
	{
	    _stats.compute(_center, _normal, this->mse, _curvature);
#ifdef DEBUG_CALC
	    _mseseq.push_back(cv::Vec2d(this->N,this->mse));
#endif
	  //nbs information to be maintained outside the class
	  //typically when initializing the graph structure
	}
#if defined(DEBUG_INIT) || defined(DEBUG_CLUSTER)
	_normalClr	= cv::Vec3b(uchar((_normal[0] + 1.0) * 0.5 * 255.0),
				    uchar((_normal[1] + 1.0) * 0.5 * 255.0),
				    uchar((_normal[2] + 1.0) * 0.5 * 255.0));
	_clr		= cv::Vec3b(rand() % 255, rand() % 255, rand() % 255);
#endif
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
	:_stats(pa._stats, pb._stats)
    {
#ifdef DEBUG_INIT
	_type = TYPE_NORMAL;
#endif
	this->nouse=false;
	this->rid = pa.N>=pb.N ? pa.rid : pb.rid;
	this->N=_stats.N;

      //ds.union(pa.rid, pb.rid) will be called later
      //in mergeNbsFrom(pa,pb) function, since
      //this object might not be accepted into the graph structure

	_stats.compute(_center, _normal, this->mse, _curvature);

#if defined(DEBUG_CLUSTER)
	_normalClr = cv::Vec3b(uchar((_normal[0] + 1.0) * 0.5 * 255.0),
			       uchar((_normal[1] + 1.0) * 0.5 * 255.0),
			       uchar((_normal[2] + 1.0) * 0.5 * 255.0));
	_clr	   = cv::Vec3b(rand() % 255, rand() % 255, rand() % 255);
#endif
      //nbs information to be maintained later if this node is accepted
    }

    void	update()
		{
		    _stats.compute(_center, _normal,
				   this->mse, _curvature);
		}

  /**
   *  \brief similarity of two plane normals
   *
   *  \param [in] p another PlaneSeg
   *  \return abs(dot(_normal, p->_normal))
   *
   *  \details 1 means identical, 0 means perpendicular
   */
    T		normalSimilarity(const PlaneSeg& p) const
		{
		    return std::abs(_normal[0] * p._normal[0] +
				    _normal[1] * p._normal[1] +
				    _normal[2] * p._normal[2]);
		}

  /**
   *  \brief signed distance between this plane and the point pt[3]
   */
    T		signedDist(const T pt[3]) const
		{
		    return _normal[0] * (pt[0] - _center[0]) +
			   _normal[1] * (pt[1] - _center[1]) +
			   _normal[2] * (pt[2] - _center[2]);
		}

  /**
   *  \brief connect this PlaneSeg to another PlaneSeg p in the graph
   *
   *  \param [in] p the other PlaneSeg
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
   *  \brief disconnect this PlaneSeg with all its neighbors
   *
   *  \details after this call, _nbs.nbs should not contain this, and _nbs should be empty i.e. after this call this PlaneSeg node should be isolated in the graph
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
   *  \param [in] pa a parent PlaneSeg of this
   *  \param [in] pb another parent PlaneSeg of this
   *  \param [in] ds the disjoint set of initial window/block membership to be updated
   *
   *  \details Only call this if this obj is accepted to be added to the graph of PlaneSeg pa and pb should not exist after this function is called, i.e. after this call this PlaneSeg node will be representing a merged node of pa and pb, and pa/pb will be isolated (and thus Garbage Collected) in the graph
   */
    void	mergeNbsFrom(PlaneSeg& pa, PlaneSeg& pb, DisjointSet& ds)
		{
		  //now we are sure that merging pa and pb is accepted
		    ds.Union(pa.rid, pb.rid);

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

		    pa.nouse=pb.nouse=true;
#ifdef DEBUG_CALC
		    if (pa.N >= pb.N)
		    {
			_mseseq.swap(pa._mseseq);
		    }
		    else
		    {
			_mseseq.swap(pb._mseseq);
		    }
		    _mseseq.push_back(cv::Vec2d(this->N,this->mse));
#endif
		}

  private:
    Stats	_stats;		//member points' 1st & 2nd order statistics
    T		_curvature;
#ifdef DEBUG_CALC
    std::vector<cv::Vec2d>	_mseseq;
#endif
#ifdef DEBUG_INIT
    enum Type
    {
	TYPE_NORMAL		= 0,	//default value
	TYPE_MISSING_DATA	= 1,
	TYPE_DEPTH_DISCONTINUE	= 2
    };

    Type	_type;
#endif
#if defined(DEBUG_INIT) || defined(DEBUG_CLUSTER)
    cv::Vec3b	_clr;
    cv::Vec3b	_normalClr;

    cv::Vec3b&	getColor(bool useNormal=true)
		{
		    if (useNormal)
			return _normalClr;
		    return _clr;
		}
#endif

  public:
    int		rid;		//root block id
    T		mse;		//mean square error
    T		_center[3]; 	//q: plane center (center of mass)
    T		_normal[3]; 	//n: plane equation n'p=q
    int		N;		//#member points, same as _stats.N
    bool	nouse;		//this PlaneSeg will be marked as nouse after merged with others to produce a new PlaneSeg node in the graph

    std::set<Ptr>	_nbs;	//neighbors, i.e. adjacency list for a graph structure
};//PlaneSeg

}//ahc
