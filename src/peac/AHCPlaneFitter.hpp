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

#include <vector>
#include <set>
#include <queue>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>

//#define EVAL_SPEED

#include "AHCTypes.hpp"
#include "AHCPlaneSeg.hpp"
#include "AHCParamSet.hpp"
#include "AHCUtils.hpp"

namespace ahc
{
using ahc::utils::Timer;
using ahc::utils::pseudocolor;


/**
 *  \brief ahc::PlaneFitter implements the Agglomerative Hierarchical Clustering based fast plane extraction
 *
 *  \details note: default parameters assume point's unit is mm
 */
template <class CLOUD>
class PlaneFitter
{
  public:
  //three types of erode operation for segmentation refinement
    enum ErodeType
    {
	ERODE_NONE	 = 0,	//no erode
	ERODE_SEG_BORDER = 1,	//erode only borders between two segments
	ERODE_ALL_BORDER = 2	//erode all borders, either between two segments or between segment and "black"
    };

  private:
    using value_t	= typename CLOUD::value_t;
    using plane_seg_t	= PlaneSeg<value_t>;
    using plane_seg_p	= typename plane_seg_t::Ptr;
    using plane_seg_sp	= typename plane_seg_t::shared_ptr;
    using param_set_t	= ParamSet<value_t>;

  //for maintaining the Min MSE heap of PlaneSeg
    struct PlaneSegMinMSECmp
    {
	bool	operator()(const plane_seg_sp& a,
			   const plane_seg_sp& b) const
		{
		    return b->_mse < a->_mse;
		}
    };

    using PlaneSegMinMSEQueue = std::priority_queue<plane_seg_sp,
						    std::vector<plane_seg_sp>,
						    PlaneSegMinMSECmp>;

  public:
  /************************************************************************/
  /* Public Class Functions                                               */
  /************************************************************************/
    PlaneFitter()
	: _cloud(nullptr),
	  _width(0),
	  _height(0),
	  _maxStep(100000),
	  _minSupport(3000),
	  _winWidth(10),
	  _winHeight(10),
	  _doRefine(true),
	  _erodeType(ERODE_ALL_BORDER),
	  _dirtyBlkMbship(true),
	  _drawCoarseBorder(false)
    {
	static const unsigned char	default_colors[10][3] =
					{
					    {255, 0, 0},
					    {255, 255, 0},
					    {100, 20, 50},
					    {0, 30, 255},
					    {10, 255, 60},
					    {80, 10, 100},
					    {0, 255, 200},
					    {10, 60, 60},
					    {255, 0, 128},
					    {60, 128, 128}
					};
	for (int i = 0; i < 10; ++i)
	{
	    _colors.push_back(cv::Vec3b(default_colors[i]));
	}
    }

    ~PlaneFitter()						{}

  /**
   *  \brief clear/reset for next run
   */
    void
    clear()
    {
	_cloud = nullptr;
	_extractedPlanes.clear();
	_ds.reset();
	_rid2plid.clear();
	_blkMap.clear();
	_rfQueue.clear();
	_dirtyBlkMbship = true;
    }

  /**
   *  \brief run AHC plane fitting on one frame of point cloud cloud
   *
   *  \param [in] _cloudIn a frame of point cloud
   *  \param [out] pMembership pointer to segmentation membership vector, each pMembership->at(i) is a vector of pixel indices that belong to the i-th extracted plane
   *  \param [out] pSeg a 3-channel RGB image as another form of output of segmentation
   *  \param [in] pIdxMap usually not needed (reserved for KinectSLAM to input pixel index map)
   *  \param [in] verbose print out cluster steps and #planes or not
   *  \return when compiled without EVAL_SPEED: 0 if _cloudIn==0 and 1 otherwise; when compiled with EVAL_SPEED: total running time for this frame
   *
   *  \details this function corresponds to Algorithm 1 in our paper
   */
    value_t
    run(const CLOUD*			cloud,
	std::vector<std::vector<int>>*	pMembership=nullptr,
	cv::Mat*			pSeg=nullptr,
	const std::vector<int>* const	pIdxMap=nullptr,
	bool				verbose=true)
    {
	if (!cloud)
	    return 0;
#ifdef EVAL_SPEED
	Timer	timer(1000), timer2(1000);
	timer.tic(); timer2.tic();
#endif
	clear();
	_cloud  = cloud;
	_height = _cloud->height();
	_width  = _cloud->width();
	_ds.reset(new DisjointSet((_height/_winHeight)*(_width/_winWidth)));

	PlaneSegMinMSEQueue	minQ;
	initGraph(minQ);
#ifdef EVAL_SPEED
	timer.toctic("init time");
#endif
	const auto	step = ahCluster(minQ);
#ifdef EVAL_SPEED
	timer.toctic("cluster time");
#endif
	if (_doRefine)
	{
	    refineDetails(pMembership, pIdxMap, pSeg);
#ifdef EVAL_SPEED
	    timer.toctic("refine time");
#endif
	}
	else
	{
	    if (pMembership)
	    {
		findMembership(*pMembership, pIdxMap);
	    }
	    if (pSeg)
	    {
		plotSegmentImage(pSeg, _minSupport);
	    }
#ifdef EVAL_SPEED
	    timer.toctic("return time");
#endif
	}
	if (verbose)
	{
	    std::cout << "#step=" << step << ", #_extractedPlanes="
		      << _extractedPlanes.size() << std::endl;
	}
#ifdef EVAL_SPEED
	return timer2.toc();
#endif
	return 1;
    }

  /**
   *  \brief print out the current parameters
   */
  // 		void logParams() const {
  // #define TMP_LOG_VAR(var) << #var "="<<(var)<<"\n"
  // 			std::cout<<"[PlaneFitter] Parameters:\n"
  // 			TMP_LOG_VAR(width)
  // 			TMP_LOG_VAR(height)
  // 			TMP_LOG_VAR(mergeMSETolerance)
  // 			TMP_LOG_VAR(initMSETolerance)
  // 			TMP_LOG_VAR(depthSigmaFactor)
  // 			TMP_LOG_VAR(similarityTh)
  //			TMP_LOG_VAR(finalMergeSimilarityTh)
  // 			TMP_LOG_VAR(simTh_znear)
  // 			TMP_LOG_VAR(simTh_zfar)
  // 			TMP_LOG_VAR(simTh_angleMin)
  // 			TMP_LOG_VAR(simTh_angleMax)
  // 			TMP_LOG_VAR(depthChangeFactor)
  // 			TMP_LOG_VAR(maxStep)
  // 			TMP_LOG_VAR(_minSupport)
  // 			TMP_LOG_VAR(_winWidth)
  // 			TMP_LOG_VAR(_winHeight)
  // 			TMP_LOG_VAR(_erodeType)
  // 			TMP_LOG_VAR(_doRefine)<<std::endl;
  // #undef TMP_LOG_VAR
  // }

  /************************************************************************/
  /* Protected Class Functions                                            */
  /************************************************************************/
  protected:
  /**
   *  \brief refine the coarse segmentation
   *
   *  \details this function corresponds to Algorithm 4 in our paper; note: plane parameters of each _extractedPlanes in the PlaneSeg is NOT updated after this call since the new points added from region grow and points removed from block erosion are not properly reflected in the PlaneSeg
   */
    void
    refineDetails(std::vector<std::vector<int>>* pMembership, //pMembership->size()==nPlanes
		  const std::vector<int>* const pIdxMap, //if pIdxMap!=0 pMembership->at(i).at(j)=pIdxMap(pixIdx)
		  cv::Mat* pSeg)
    {
	if (!pMembership && !pSeg)
	    return;

	std::vector<bool>	isValidExtractedPlane; //some planes might be eroded completely
	findBlockMembership(isValidExtractedPlane);

      //save the border regions
	std::vector<int>	border;
	if (_drawCoarseBorder && pSeg)
	{
	    border.resize(_rfQueue.size());
	    for (int i = 0; i < (int)_rfQueue.size(); ++i)
	    {
		border[i] = _rfQueue[i].first;
	    }
	}

	floodFill();

      //try to merge one last time
	std::vector<plane_seg_sp>	oldExtractedPlanes;
	_extractedPlanes.swap(oldExtractedPlanes);
	PlaneSegMinMSEQueue		minQ;
	for (int i = 0; i < (int)oldExtractedPlanes.size(); ++i)
	{
	    if (isValidExtractedPlane[i])
		minQ.push(oldExtractedPlanes[i]);
	}
	ahCluster(minQ, false);

      //find plane idx maping from oldExtractedPlanes to final _extractedPlanes
	std::vector<int>	plidmap(oldExtractedPlanes.size(), -1);
	int			nFinalPlanes = 0;
	for (int i = 0; i < (int)oldExtractedPlanes.size(); ++i)
	{
	    if (!isValidExtractedPlane[i])
	    {
		plidmap[i] = -1;//this plane was eroded
		continue;
	    }

	    const auto&	op     = *oldExtractedPlanes[i];
	    const auto	np_rid = _ds->Find(op._rid);
	    if (np_rid == op._rid)
	    {//op is not changed
		if (plidmap[i] < 0) //if not updated, otherwise was updated before
		    plidmap[i] = nFinalPlanes++;
	    }
	    else
	    {//op is merged to np
		const auto	npid = _rid2plid[np_rid];
		if (plidmap[npid] < 0) //if np's idmap not updated
		    plidmap[i] = plidmap[npid] = nFinalPlanes++;
		else
		    plidmap[i] = plidmap[npid];
	    }
	}
	assert(nFinalPlanes == (int)_extractedPlanes.size());

      //scan _membershipImg
	if (nFinalPlanes > _colors.size())
	{
	    auto	tmpColors = pseudocolor(nFinalPlanes -
						(int)_colors.size());
	    _colors.insert(_colors.end(), tmpColors.begin(), tmpColors.end());
	}
	if (pMembership)
	{
	    pMembership->resize(nFinalPlanes, std::vector<int>());
	    for (int i = 0; i < nFinalPlanes; ++i)
	    {
		pMembership->at(i).reserve(
		    (int)(_extractedPlanes[i]->_N * 1.2f));
	    }
	}

	static const cv::Vec3b	blackColor(0, 0, 0);
	const auto		nPixels = _width*_height;
	for (int i = 0; i < nPixels; ++i)
	{
	    auto&	plid=_membershipImg.at<int>(i);
	    if (plid >= 0 && plidmap[plid] >= 0)
	    {
		plid = plidmap[plid];
		if (pSeg)
		    pSeg->at<cv::Vec3b>(i) = _colors[plid];
		if (pMembership)
		    pMembership->at(plid).push_back(pIdxMap ? pIdxMap->at(i)
							    : i);
	    }
	    else
	    {
		if (pSeg)
		    pSeg->at<cv::Vec3b>(i) = blackColor;
	    }
	}

	static const cv::Vec3b	whiteColor(255, 255, 255);
	for (int k = 0; pSeg && _drawCoarseBorder && k < (int)border.size();
	     ++k)
	{
	    pSeg->at<cv::Vec3b>(border[k]) = whiteColor;
	}
      //TODO: refine the plane equation as well after!!
    }

  /**
   *  \brief find out all valid 4-connect neighbours pixels of pixel (i,j)
   *
   *  \param [in] i row index of the center pixel
   *  \param [in] j column index of the center pixel
   *  \param [in] H height of the image
   *  \param [in] W weight of the image
   *  \param [out] nbs pixel id of all valid neighbours
   *  \return number of valid neighbours
   *
   *  \details invalid 4-connect neighbours means out of image boundary
   */
    static int
    getValid4Neighbor(const int i, const int j,
		      const int H, const int W, int nbs[4])
    {
	const auto	id  = i*W + j;
	int		cnt = 0;
	if (j > 0)
	    nbs[cnt++] = id - 1;		//left
	if (j < W - 1)
	    nbs[cnt++] = id + 1;		//right
	if (i > 0)
	    nbs[cnt++] = id - W;		//up
	if (i < H - 1)
	    nbs[cnt++] = id + W;		//down
	return cnt;
    }

  /**
   *  \brief find out pixel (pixX, pixY) belongs to which initial block/window
   *
   *  \param [in] pixX column index
   *  \param [in] pixY row index
   *  \return initial block id, or -1 if not in any block (usually because _winWidth%_width!=0 or _winHeight%height!=0)
   */
    inline int
    getBlockIdx(const int pixX, const int pixY) const
    {
	assert(pixX >= 0 && pixY >= 0 && pixX < _width && pixY < _height);
	const auto	Nw = _width/_winWidth;
	const auto	Nh = _height/_winHeight;
	const auto	by = pixY/_winHeight;
	const auto	bx = pixX/_winWidth;
	return (by < Nh && bx < Nw) ? (by*Nw + bx) : -1;
    }

  /**
   *  \brief region grow from coarse segmentation boundaries
   *
   *  \details this function implemented line 14~25 of Algorithm 4 in our paper
   */
    void
    floodFill()
    {
	std::vector<float> distMap(_height*_width,
				   std::numeric_limits<float>::max());

	for (int k=0; k<(int)_rfQueue.size(); ++k)
	{
	    const auto	sIdx  = _rfQueue[k].first;
	    const auto	seedy = sIdx/_width;
	    const auto	seedx = sIdx-seedy*_width;
	    const auto	plid  =_rfQueue[k].second;
	    const auto&	pl    = *_extractedPlanes[plid];

	    int		nbs[4] = {-1};
	    const auto	Nnbs=getValid4Neighbor(seedy, seedx,
					       _height, _width, nbs);
	    for (int itr = 0; itr < Nnbs; ++itr)
	    {
		const int	cIdx  = nbs[itr];
		int&		trail = _membershipImg.at<int>(cIdx);
		if (trail<=-6)
		    continue; //visited from 4 neighbors already, skip
		if (trail>=0 && trail==plid)
		    continue; //if visited by the same plane, skip

		const auto	cy    = cIdx/_width;
		const auto	cx    = cIdx-cy*_width;
		const auto	blkid = getBlockIdx(cx, cy);
		if (blkid >= 0 && _blkMap[blkid] >= 0)
		    continue; //not in "black" block

		value_t	pt[3] = {0};
		float	cdist = -1;
		if (_cloud->get(cy, cx, pt[0], pt[1], pt[2]) &&
		    std::pow(cdist = float(std::abs(pl.signedDist(pt))), 2)
		    < 9*pl._mse + 1e-5) //point-plane distance within 3*std
		{
		    if (trail >= 0)
		    {
			auto&	n_pl = *_extractedPlanes[trail];
			if (pl.normalSimilarity(n_pl) >=
			    _params.T_ang(param_set_t::P_REFINE,
					  pl._center[2]))
			{//potential for merging
			    n_pl.connect(_extractedPlanes[plid].get());
			}
		    }
		    auto&	old_dist = distMap[cIdx];
		    if (cdist < old_dist)
		    {
			trail	 = plid;
			old_dist = cdist;
			_rfQueue.push_back(std::make_pair(cIdx, plid));
		    }
		    else if (trail < 0)
		    {
			trail -= 1;
		    }
		}
		else
		{
		    if (trail<0)
			trail -= 1;
		}
	    }
	}//for _rfQueue
    }

  /**
   *  \brief erode each segment at initial block/window level
   *
   *  \param [in] isValidExtractedPlane coarsely extracted plane i is completely eroded if isValidExtractedPlane(i)==false
   *
   *  \details this function implements line 5~13 of Algorithm 4 in our paper, called by refineDetails; FIXME: after this ds is not updated, i.e. is dirtied
   */
    void
    findBlockMembership(std::vector<bool>& isValidExtractedPlane)
    {
	_rid2plid.clear();
	for (int plid = 0; plid < (int)_extractedPlanes.size(); ++plid)
	{
	    _rid2plid.insert(std::make_pair(_extractedPlanes[plid]->_rid,
					    plid));
	}

	const auto	Nh	   = _height/_winHeight;
	const auto	Nw	   = _width/_winWidth;
	const auto	NptsPerBlk = _winHeight*_winWidth;

	_membershipImg.create(_height, _width, CV_32SC1);
	_membershipImg.setTo(-1);
	_blkMap.resize(Nh*Nw);

	isValidExtractedPlane.resize(_extractedPlanes.size(), false);
	for (int i = 0,blkid = 0; i < Nh; ++i)
	{
	    for (int j = 0; j < Nw; ++j, ++blkid)
	    {
		const auto	setid   = _ds->Find(blkid);
		const auto	setSize = _ds->getSetSize(setid)*NptsPerBlk;
		if (setSize>=_minSupport)
		{//cluster large enough
		    int		nbs[4] = {-1};
		    const auto	nNbs = getValid4Neighbor(i, j, Nh, Nw, nbs);
		    bool	nbClsAllTheSame = true;
		    for (int k = 0; k < nNbs && _erodeType != ERODE_NONE; ++k)
		    {
			if (_ds->Find(nbs[k]) != setid &&
			    (_erodeType == ERODE_ALL_BORDER ||
			     _ds->getSetSize(nbs[k])*NptsPerBlk >=
			     _minSupport))
			{
			    nbClsAllTheSame = false;
			    break;
			}
		    }
		    const auto	plid = _rid2plid[setid];
		    if (nbClsAllTheSame)
		    {
			_blkMap[blkid]=plid;
			const auto	by = blkid/Nw;
			const auto	bx = blkid - by*Nw;
			_membershipImg(
			    cv::Range(by*_winHeight, (by + 1)*_winHeight),
			    cv::Range(bx*_winWidth,  (bx + 1)*_winWidth))
			    .setTo(plid);
			isValidExtractedPlane[plid] = true;
		    }
		    else
		    {//erode border region
			_blkMap[blkid] = -1;
		      //_extractedPlanes[plid]->stats.pop(blkStats[blkid]);
		    }
		}
		else
		{//too small cluster, i.e. "black" cluster
		    _blkMap[blkid] = -1;
		}//if setSize>=blkMinSupport

	      //save seed points for floodFill
		if (_blkMap[blkid] <0 )
		{//current block is not valid
		    if (i > 0)
		    {
			const auto	u_blkid = blkid - Nw;
			if (_blkMap[u_blkid] >= 0)
			{//up blk is in border
			    const auto	u_plid  =_blkMap[u_blkid];
			    const auto	spixidx = (i*_winHeight - 1)*_width
						+ j*_winWidth;
			    for (int k = 1; k < _winWidth; ++k)
			    {
				_rfQueue.push_back(
				    std::make_pair(spixidx + k, u_plid));
			    }
			}
		    }
		    if (j > 0)
		    {
			const auto	l_blkid = blkid-1;
			if (_blkMap[l_blkid] >= 0)
			{//left blk is in border
			    const auto	l_plid  =_blkMap[l_blkid];
			    const auto	spixidx = (i*_winHeight)*_width
						+ j*_winWidth - 1;
			    for (int k = 0; k < _winHeight - 1; ++k)
			    {
				_rfQueue.push_back(
				    std::make_pair(spixidx + k*_width, l_plid));
			    }
			}
		    }
		}
		else
		{//current block is still valid
		    const int plid=_blkMap[blkid];
		    if (i>0)
		    {
			const auto	u_blkid = blkid-Nw;
			if (_blkMap[u_blkid] != plid)
			{//up blk is in border
			    const auto	spixidx = (i*_winHeight)*_width
						+ j*_winWidth;
			    for (int k = 0; k < _winWidth - 1; ++k)
			    {
				_rfQueue.push_back(
				    std::make_pair(spixidx +k, plid));
			    }
			}
		    }
		    if (j > 0)
		    {
			const auto	l_blkid = blkid - 1;
			if (_blkMap[l_blkid]!=plid)
			{//left blk is in border
			    const auto	spixidx = (i*_winHeight)*_width
						+ j*_winWidth;
			    for (int k = 1; k < _winHeight; ++k)
			    {
				_rfQueue.push_back(
				    std::make_pair(spixidx + k*_width, plid));
			    }
			}
		    }
		}//save seed points for floodFill
	    }
	}//for blkik

      ////update plane equation
      //for (int i=0; i<(int)_extractedPlanes.size(); ++i) {
      //	if (isValidExtractedPlane[i]) {
      //		if (_extractedPlanes[i]->stats.N>=_minSupport)
      //			_extractedPlanes[i]->update();
      //	} else {
      //		_extractedPlanes[i]->_nouse=true;
      //	}
      //}
    }

  //called by findMembership and/or plotSegmentImage when _doRefine==false
    void
    findBlockMembership()
    {
	if (!_dirtyBlkMbship)
	    return;
	_dirtyBlkMbship = false;

	_rid2plid.clear();
	for (int plid=0; plid < (int)_extractedPlanes.size(); ++plid)
	{
	    _rid2plid.insert(std::make_pair(_extractedPlanes[plid]->_rid,
					    plid));
	}

	const auto Nh	      = _height/_winHeight;
	const auto Nw	      = _width/_winWidth;
	const auto NptsPerBlk = _winHeight*_winWidth;

	_membershipImg.create(_height, _width, CV_32SC1);
	_membershipImg.setTo(-1);
	_blkMap.resize(Nh*Nw);

	for (int i = 0, blkid = 0; i < Nh; ++i)
	{
	    for (int j = 0; j < Nw; ++j, ++blkid)
	    {
		const auto	setid   = _ds->Find(blkid);
		const auto	setSize = _ds->getSetSize(setid)*NptsPerBlk;
		if (setSize >= _minSupport)
		{//cluster large enough
		    const auto	plid=_rid2plid[setid];
		    _blkMap[blkid]=plid;
		    const auto	by = blkid/Nw;
		    const auto	bx = blkid - by*Nw;
		    _membershipImg(
			cv::Range(by*_winHeight, (by + 1)*_winHeight),
			cv::Range(bx*_winWidth,  (bx + 1)*_winWidth))
			.setTo(plid);
		}
		else
		{//too small cluster, i.e. "black" cluster
		    _blkMap[blkid]=-1;
		}//if setSize>=blkMinSupport
	    }
	}//for blkik
    }

		//called by run when _doRefine==false
    void
    findMembership(std::vector< std::vector<int> >& ret,
		   const std::vector<int>* pIdxMap)
    {
	const auto	Nh = _height/_winHeight;
	const auto	Nw = _width/_winWidth;
	findBlockMembership();
	const auto	cnt = (int)_extractedPlanes.size();
	ret.resize(cnt);
	for (int i = 0; i < cnt; ++i)
	    ret[i].reserve(_extractedPlanes[i]->_N);
	for (int i = 0, blkid = 0; i < Nh; ++i)
	{
	    for (int j = 0; j < Nw; ++j, ++blkid)
	    {
		const auto	plid = _blkMap[blkid];
		if (plid < 0)
		    continue;
		for (int y = i*_winHeight; y < (i + 1)*_winHeight; ++y)
		{
		    for (int x = j*_winWidth; x < (j + 1)*_winWidth; ++x)
		    {
			const auto	pixIdx = x + y*_width;
			ret[plid].push_back(pIdxMap ? pIdxMap->at(pixIdx)
						    : pixIdx);
		    }
		}
	    }
	}
    }

		//called by run when _doRefine==false
    void
    plotSegmentImage(cv::Mat* pSeg, const value_t supportTh)
    {
	if (pSeg==0)
	    return;

	const auto	Nh  = _height/_winHeight;
	const auto	Nw  = _width/_winWidth;
	int		cnt = 0;

	std::vector<int>	ret;
	std::vector<int>*	pBlkid2plid;
	if (supportTh==_minSupport)
	{
	    findBlockMembership();
	    pBlkid2plid = &(_blkMap);
	    cnt		= (int)_extractedPlanes.size();
	}
	else
	{ //mainly for DEBUG_CLUSTER since then supportTh!=_minSupport
	    std::map<int, int>	map; //map setid->cnt
	    ret.resize(Nh*Nw);

	    for (int i = 0, blkid = 0; i < Nh; ++i)
	    {
		for (int j = 0; j < Nw; ++j, ++blkid)
		{
		    const auto	setid   = _ds->Find(blkid);
		    const auto	setSize = _ds->getSetSize(setid)
					* _winHeight*_winWidth;
		    if (setSize >= supportTh)
		    {
			auto	fitr = map.find(setid);
			if (fitr == map.end())
			{//found a new set id
			    map.insert(std::make_pair(setid, cnt));
			    ret[blkid] = cnt;
			    ++cnt;
			}
			else
			{//found a existing set id
			    ret[blkid] = fitr->second;
			}
		    }
		    else
		    {//too small cluster, ignore
			ret[blkid] = -1;
		    }
		}
	    }

	    pBlkid2plid = &ret;
	}
	std::vector<int>&	blkid2plid = *pBlkid2plid;

	if (cnt > _colors.size())
	{
	    auto	tmpColors = pseudocolor(cnt - (int)_colors.size());
	    _colors.insert(_colors.end(), tmpColors.begin(), tmpColors.end());
	}
	auto&	seg = *pSeg;
	for (int i = 0, blkid = 0; i < Nh; ++i)
	{
	    for (int j = 0; j < Nw; ++j, ++blkid)
	    {
		const auto	plid = blkid2plid[blkid];
		if (plid >= 0)
		{
		    seg(cv::Range(i*_winHeight, (i + 1)*_winHeight),
			cv::Range(j*_winWidth,  (j + 1)*_winWidth))
			.setTo(_colors[plid]);
		}
		else
		{
		    seg(cv::Range(i*_winHeight, (i + 1)*_winHeight),
			cv::Range(j*_winWidth,  (j + 1)*_winWidth))
			.setTo(cv::Vec3b(0,0,0));
		}
	    }
	}
    }

  /**
   *  \brief initialize a graph from cloud
   *
   *  \param [in/out] minQ a min MSE queue of PlaneSegs
   *
   *  \details this function implements Algorithm 2 in our paper
   */
    void
    initGraph(PlaneSegMinMSEQueue& minQ)
    {
	const auto	Nh = _height/_winHeight;
	const auto	Nw = _width/_winWidth;

      //1. init nodes
	std::vector<plane_seg_p>	G(Nh*Nw, 0);
      //blkStats.resize(Nh*Nw);

	for (int i = 0; i < Nh; ++i)
	{
	    for (int j = 0; j < Nw; ++j)
	    {
		plane_seg_sp	p(new plane_seg_t(*_cloud,
						  i*Nw + j,
						  i*_winHeight, j*_winWidth,
						  _width,	_height,
						  _winWidth,    _winHeight,
						  _params));
		if (p->_mse<_params.T_mse(param_set_t::P_INIT, p->_center[2])
		    && !p->_nouse)
		{
		    G[i*Nw+j] = p.get();
		    minQ.push(p);
		  //blkStats[i*Nw+j]=p->stats;
		}
		else
		{
		    G[i*Nw+j] = 0;
		}
	    }
	}

      //2. init edges
      //first pass, connect neighbors from row direction
	for (int i = 0; i < Nh; ++i)
	{
	    for (int j = 1; j < Nw; j+=2)
	    {
		const auto	cidx = i*Nw + j;
		if (G[cidx - 1] == 0)
		{
		    --j;
		    continue;
		}
		if (G[cidx] == 0)
		    continue;
		if (j < Nw-1 && G[cidx+1] == 0)
		{
		    ++j;
		    continue;
		}

		const auto similarityTh = _params.T_ang(param_set_t::P_INIT,
							G[cidx]->_center[2]);
		if ((j < Nw - 1 &&
		     G[cidx - 1]->normalSimilarity(*G[cidx + 1]) >=
		     similarityTh) ||
		   (j == Nw - 1 &&
		    G[cidx]->normalSimilarity(*G[cidx - 1]) >= similarityTh))
		{
		    G[cidx]->connect(G[cidx - 1]);
		    if (j < Nw - 1)
			G[cidx]->connect(G[cidx + 1]);
		}
		else
		{//otherwise current block is in edge region
		    --j;
		}
	    }
	}
      //second pass, connect neighbors from column direction
	for (int j = 0; j < Nw; ++j)
	{
	    for (int i = 1; i < Nh; i += 2)
	    {
		const auto	cidx = i*Nw + j;
		if (G[cidx-Nw] == 0)
		{
		    --i;
		    continue;
		}
		if (G[cidx] == 0)
		    continue;
		if (i < Nh - 1 && G[cidx+Nw] == 0)
		{
		    ++i;
		    continue;
		}

		const auto similarityTh=_params.T_ang(param_set_t::P_INIT,
						      G[cidx]->_center[2]);
		if ((i < Nh - 1 &&
		     G[cidx - Nw]->normalSimilarity(*G[cidx+Nw]) >=
		     similarityTh) ||
		    (i == Nh - 1 &&
		     G[cidx]->normalSimilarity(*G[cidx-Nw]) >= similarityTh))
		{
		    G[cidx]->connect(G[cidx - Nw]);
		    if (i < Nh - 1)
			G[cidx]->connect(G[cidx + Nw]);
		}
		else
		{
		    --i;
		}
	    }
	}
    }

  /**
   *  \brief main clustering step
   *
   *  \param [in] minQ a min MSE queue of PlaneSegs
   *  \param [in] debug whether to collect some statistics when compiled with DEBUG_CALC
   *  \return number of cluster steps
   *
   *  \details this function implements the Algorithm 3 in our paper
   */
    int
    ahCluster(PlaneSegMinMSEQueue& minQ, bool debug=true)
    {
	int step = 0;
	while (!minQ.empty() && step <= _maxStep)
	{
	    plane_seg_sp	p = minQ.top();
	    minQ.pop();
	    if (p->_nouse)
	    {
		assert(p->_nbs.size() <= 0);
		continue;
	    }

	    plane_seg_sp	cand_merge;
	    plane_seg_p		cand_nb(0);
	    for (auto nb : p->_nbs)
	    {//test merge with all nbs, pick the one with min mse
	      //TODO: should we use dynamic similarityTh here?
	      //const value_t similarityTh=ahc::depthDependNormalDeviationTh(p->_center[2],500,4000,M_PI*15/180.0,M_PI/2);
		if (p->normalSimilarity(*nb) <
		    _params.T_ang(param_set_t::P_MERGING, p->_center[2]))
		    continue;
		plane_seg_sp	merge(new plane_seg_t(*p, *nb));
		if (cand_merge == 0 || cand_merge->_mse > merge->_mse ||
		   (cand_merge->_mse == merge->_mse &&
		    cand_merge->_N < merge->_mse))
		{
		    cand_merge=merge;
		    cand_nb=nb;
		}
	    }//for nbs

	  //TODO: maybe a better merge condition? such as adaptive threshold on MSE like Falzenszwalb's method
	    if (cand_merge!=0 &&
		cand_merge->_mse < _params.T_mse(param_set_t::P_MERGING,
						 cand_merge->_center[2]))
	    {//merge and add back to minQ
		minQ.push(cand_merge);
		cand_merge->mergeNbsFrom(*p, *cand_nb, *_ds);
	    }
	    else
	    {//do not merge, but extract p
		if (p->_N >= _minSupport)
		{
		    _extractedPlanes.push_back(p);
		}
		p->disconnectAllNbs();
	    }
	    ++step;
	}//end while minQ

	while (!minQ.empty())
	{//just check if any remaining PlaneSeg if maxstep reached
	    const auto	p = minQ.top();
	    minQ.pop();
	    if (p->_N >= _minSupport)
	    {
		_extractedPlanes.push_back(p);
	    }
	    p->disconnectAllNbs();
	}
      //static PlaneSegSizeCmp sizecmp;
	std::sort(_extractedPlanes.begin(), _extractedPlanes.end(),
		  [](const auto& a, const auto& b){ return b->_N < a->_N; });
	return step;
    }

  private:
  // [Input]
  //! dim=<heightxwidthx3>, no ownership
    const CLOUD*			_cloud;
  //! witdth=#cols, _height=#rows (size of the input point cloud)
    int					_width, _height;
  //! max number of steps for merging clusters
    int					_maxStep;
  //! min number of supporting point
    int					_minSupport;
  //! make sure _width is divisible by _winWidth
    int					_winWidth;
  //! similarly for _height and _winHeight
    int					_winHeight;
  //! perform refinement of details or not
    bool				_doRefine;
    ErodeType				_erodeType;
  //! sets of parameters controlling dynamic thresholds T_mse, T_ang, T_dz
    param_set_t				_params;

  // [Output]
  //! with ownership, this disjoint set maintains membership of initial window/blocks during AHC merging
    ahc::shared_ptr<DisjointSet>	_ds;
  //! a set of extracted planes
    std::vector<plane_seg_sp>		_extractedPlanes;
  //! segmentation map of the input pointcloud, membershipImg(i,j) records which plane (plid, i.e. plane id) this pixel/point (i,j) belongs to
    cv::Mat				_membershipImg;

  // [Intermediate]
  //! _extractedPlanes[_rid2plid[rootid]].rid==rootid, i.e. _rid2plid[rid] gives the idx of a plane in _extractedPlanes
    std::map<int, int>			_rid2plid;
  //! (i,j) block belong to _extractedPlanes[blkMap[i*Nh+j]]
    std::vector<int>			_blkMap;
  //! std::vector<std::vector<int>> blkMembership; //blkMembership[i] contains all block id for _extractedPlanes[i]
    bool				_dirtyBlkMbship;
    std::vector<cv::Vec3b>		_colors;
  //for region grow/floodfill, p.first=pixidx, p.second=plid
    std::vector<std::pair<int,int>>	_rfQueue;
    bool				_drawCoarseBorder;
  //std::vector<plane_seg_t::Stats> blkStats;
};//end of PlaneFitter
}//end of namespace ahc
