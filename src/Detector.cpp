/*!
* \file		Detector.cpp
* \author	Toshio UESHIBA
*/
#include "Detector.h"
#include <aist_utility/sensor_msgs.h>
#include <aist_utility/opencv.h>

namespace plane_detector
{
/************************************************************************
*  global functions							*
************************************************************************/
template <class T> inline T
val(const sensor_msgs::Image& image_msg, int u, int v)
{
    using namespace	sensor_msgs;

    if (image_msg.encoding == image_encodings::TYPE_16UC1)
    	return T(0.001) * *reinterpret_cast<const uint16_t*>(
    				image_msg.data.data() + v*image_msg.step
    						      + u*sizeof(uint16_t));
    else
	return *reinterpret_cast<const T*>(image_msg.data.data()
					   + v*image_msg.step + u*sizeof(T));
}

/************************************************************************
*  class Detector							*
************************************************************************/
Detector::Detector(const ros::NodeHandle& nh)
    :_nh(nh),
     _it(_nh),
     _camera_info_sub(_nh, "/camera_info", 1),
     _image_sub(_it, "/image", 1),
     _depth_sub(_it, "/depth", 1),
     _sync(sync_policy_t(10), _camera_info_sub, _image_sub, _depth_sub),
     _image_pub(_it.advertise("result", 1)),
     _pose_pub(_nh.advertise<geometry_msgs::PoseStamped>("pose", 100)),
     _ddr(),
     _planarityTolerance(0.001),
     _cloud(),
     _plane_fitter(),
     _plane_vertices()
{
  // Set planarity tolerance and setup its ddynamic_recoconfigure service.
    _ddr.registerVariable<double>(
    	"planarity_tolerance", &_planarityTolerance,
    	"Planarity tolerance for extracting marker region(in meters)",
    	0.0005, 0.05);

  // Pulish ddynamic_reconfigure service.
    _ddr.publishServicesTopics();

  // Register callback for marker detection.
    _sync.registerCallback(&Detector::detect_plane_cb, this);
}

void
Detector::run()
{
    ros::spin();
}

void
Detector::detect_plane_cb(const camera_info_p& camera_info_msg,
			  const image_p& image_msg, const image_p& depth_msg)
{
    try
    {
	using namespace	sensor_msgs;
	using namespace aist_utility;

	_cloud.resize(depth_msg->height, depth_msg->width);
	depth_to_points<float>(*camera_info_msg, *depth_msg, _cloud.begin(),
			       milimeters<float>);

	_seg_img.header	  = depth_msg->header;
	_seg_img.encoding = image_encodings::RGB8;
	if (_seg_img.image.rows != depth_msg->height ||
	    _seg_img.image.cols != depth_msg->width)
	    _seg_img.image = cv::Mat(depth_msg->height, depth_msg->width,
				     CV_8UC3);

	_plane_fitter.run(&_cloud, &_plane_vertices, &_seg_img.image);

	_image_pub.publish(_seg_img.toImageMsg());
    }
    catch (const std::exception& e)
    {
	ROS_WARN_STREAM(e.what());
    }
}

template <class T> inline cv::Vec<T, 3>
Detector::at(const image_t& depth_msg, int u, int v) const
{
    const auto	xyz = view_vector<T>(u, v);
    const auto	d   = val<T>(depth_msg, u, v);

    return {xyz[0]*d, xyz[1]*d, d};
}

template <class T> cv::Vec<T, 3>
Detector::at(const image_t& depth_msg, T u, T v) const
{
    const int	u0 = std::floor(u);
    const int	v0 = std::floor(v);
    const int	u1 = std::ceil(u);
    const int	v1 = std::ceil(v);
    for (auto vv = v0; vv <= v1; ++vv)
	for (auto uu = u0; uu <= u1; ++uu)
	{
	    const auto	xyz = at<T>(depth_msg, uu, vv);
	    if (xyz[2] != T(0))
		return xyz;
	}

    return {T(0), T(0), T(0)};
}

}	// namespace plane_detector
