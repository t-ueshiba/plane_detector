/*
 *  \file	Parameters.h
 *  \author	Toshio UESHIBA
 */
#pragma once

#include <cmath>

namespace plane_detector
{
template <class T>
struct Parameters
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
    Parameters(InitType initType=INIT_STRICT)
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

}	// namespace plane_detector
