/*!
 *  \file	nodelet.cpp
 *  \author	Toshio UESHIBA
 */
#include "Detector.h"
#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>

namespace plane_detector
{
/************************************************************************
*  class DetectorNodelet						*
************************************************************************/
class DetectorNodelet : public nodelet::Nodelet
{
  public:
			DetectorNodelet()				{}

    virtual void	onInit()					;

  private:
    boost::shared_ptr<Detector>	_node;
};

void
DetectorNodelet::onInit()
{
    NODELET_INFO("plane_detector::DetectorNodelet::onInit()");
    _node.reset(new Detector(getPrivateNodeHandle()));
}

}	// namespace plane_detector

PLUGINLIB_EXPORT_CLASS(plane_detector::DetectorNodelet, nodelet::Nodelet);
