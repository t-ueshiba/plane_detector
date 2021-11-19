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
     _camera_info_pub(_nh.advertise<camera_info_t>("camera_info", 1)),
     _image_pub(_it.advertise("image", 1)),
     _depth_pub(_it.advertise("depth", 1)),
     _pose_pub(_nh.advertise<geometry_msgs::PoseStamped>("pose", 100)),
     _cloud_sub(_nh.subscribe<cloud_t>("/pointcloud", 1,
				       boost::bind(&Detector::cloud_cb,
						   this, _1))),
     _cloud(),
     _plane_fitter(),
     _plane_vertices()
{
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

	_camera_info_pub.publish(camera_info_msg);
	_image_pub.publish(_seg_img.toImageMsg());
	_depth_pub.publish(depth_msg);
    }
    catch (const std::exception& e)
    {
	ROS_WARN_STREAM(e.what());
    }
}

void
Detector::cloud_cb(const cloud_p& cloud_msg)
{
    try
    {
	using namespace	sensor_msgs;
	using namespace	aist_utility;
	
	_cloud.resize(cloud_msg->height, cloud_msg->width);
	cloud_to_points<float>(*cloud_msg, _cloud.begin(), milimeters<float>);
	
	_seg_img.header	  = cloud_msg->header;
	_seg_img.encoding = image_encodings::RGB8;
	if (_seg_img.image.rows != cloud_msg->height ||
	    _seg_img.image.cols != cloud_msg->width)
	    _seg_img.image = cv::Mat(cloud_msg->height, cloud_msg->width,
				     CV_8UC3);

	_plane_fitter.run(&_cloud, &_plane_vertices, &_seg_img.image);

	_image_pub.publish(_seg_img.toImageMsg());
    }
    catch (const std::exception& e)
    {
	ROS_WARN_STREAM(e.what());
    }
}
    
}	// namespace plane_detector
