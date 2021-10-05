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

#include <cmath>

namespace ahc
{
#define MACRO_DEG2RAD(d) ((d)*M_PI/180.0)
#define MACRO_RAD2DEG(r) ((r)*180.0/M_PI)

/**
*  \brief	ParamSet is a struct representing a set of parameters used
*		in ahc::PlaneFitter
*/
template <class T>
struct ParamSet
{
  public:
    enum InitType
    {
	INIT_STRICT = 0,  // no nan point is allowed in any valid init blocks
	INIT_LOOSE  = 1   // at most half of a init block can be nan point
    };

    enum Phase
    {
	P_INIT	  = 0,
	P_MERGING = 1,
	P_REFINE  = 2
    };

  public:
    ParamSet(InitType initType=INIT_STRICT)
	:_initType(initType),
	 _depthSigma(1.6e-6),
	 _stdTol_init(5),
	 _stdTol_merge(8),
	 _z_near(500),
	 _z_far(4000),
	 _angle_near(MACRO_DEG2RAD(15.0)),
	 _angle_far(MACRO_DEG2RAD(90.0)),
	 _similarityTh_merge(std::cos(MACRO_DEG2RAD(60.0))),
	 _similarityTh_refine(std::cos(MACRO_DEG2RAD(30.0))),
	 _depthAlpha(0.04),
	 _depthChangeTol(0.02)
    {
    }

    InitType	initType()			const	{ return _initType; }

  /**
   *  \brief	Dynamic MSE threshold, depending on depth z
   *
   *  \param	[in] phase specify whether invoked in initGraph or when merging
   *  \param	[in] z current depth, unit: u mm
   *  \return	the MSE threshold at depth z, unit: u^2 mm^2
   *
   *  \details	Reference: 2012.Sensors.Khoshelham.Accuracy and Resolution
   *		of Kinect Depth Data for Indoor Mapping Applications
   */
    T	T_mse(Phase phase, T z=0) const
	{
	  // theoretical point-plane distance std = sigma * z * z
	  // sigma corresponds to \sigma * (m/f/b) in the 2012.Khoshelham paper
	  // we add a stdTol to move the theoretical curve up as tolerances
	    switch (phase)
	    {
	      case P_INIT:
		return std::pow(_depthSigma * z * z + _stdTol_init, 2);
	      case P_MERGING:
	      case P_REFINE:
	      default:
		break;
	    }

	    return std::pow(_depthSigma * z * z + _stdTol_merge, 2);
	}

  /**
   *  \brief	Dynamic normal deviation threshold, depending on depth z
   *
   *  \param	[in] phase specify whether invoked in initGraph or when merging
   *  \param	[in] z current depth (z>=0)
   *  \return	cos of the normal deviation threshold at depth z
   *
   *  \details	This is simply a linear mapping from depth to thresholding
   *		angle and the threshold will be used to reject edge
   *		when initialize the graph;
   *		this function corresponds to T_{ANG} in our paper
   */
    T	T_ang(Phase phase, T z=0) const
	{
	    switch (phase)
	    {
	      case P_INIT:
	      {//linear maping z->thresholding angle, clipping z also
		z = std::min(std::max(z, _z_near), _z_far);
		const T	factor = (_angle_far - _angle_near)/(_z_far - _z_near);
		return std::cos(factor * (z- _z_near) + _angle_near);
	      }
	      case P_MERGING:
		return _similarityTh_merge;
	      default:
		break;
	    }

	    return _similarityTh_refine;
	}

  /**
   *  \brief	Dynamic threshold to test whether the two adjacent pixels are
   *		discontinuous in depth
   *
   *  \param	[in] z depth of the current pixel
   *  \return	the max depth change allowed at depth z
   *		to make the points connected in a single block
   *
   *  \details	This is modified from pcl's segmentation code as well as
   *		suggested in 2013.iros.Holzer essentially returns
   *		factor * z + tolerance
   *		(TODO: maybe change this to 3D-point distance threshold)
   */
    T	T_dz(T z) const
	{
	    return _depthAlpha * std::abs(z) + _depthChangeTol;
	}

  private:
    const InitType	_initType;

  /* related to T_mse	*/
  //! \sigma in the paper, unit: u^-1 mm^-1
    const T		_depthSigma;
  //! \epsilon in the paper, used when init graph, unit: u mm
    const T		_stdTol_init;
  //! \epsilon in the paper, used when merging nodes, unit: u mm
    const T		_stdTol_merge;

  /* related to T_ang	*/
  //! unit: u mm, closest/farthest z to be considered
    const T		_z_near, _z_far;
  //! unit: rad, corresponding normal deviation angle
    const T		_angle_near, _angle_far;
  //! unit: none, 1 means the same, 0 means perpendicular
    const T		_similarityTh_merge;
    const T		_similarityTh_refine;	//!< unit: none

  /* related to T_dz	*/
  //! unit: none, corresponds to the 2*\alpha in the paper
    const T		_depthAlpha;
    const T		_depthChangeTol;	//!< unit: u mm
};

}	// namespace ahc
