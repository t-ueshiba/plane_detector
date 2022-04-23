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
*  class Detector							*
************************************************************************/
Detector::Detector(const ros::NodeHandle& nh)
    :_nh(nh),
     _it(_nh),
     _camera_sub(_it.subscribeCamera("/depth", 1, &Detector::camera_cb, this)),
     _cloud_sub(_nh.subscribe<cloud_t>("/pointcloud", 1,
				       boost::bind(&Detector::cloud_cb,
						   this, _1))),
     _image_pub(_it.advertise("image", 1)),
     _cloud(),
     _plane_fitter(),
     _plane_vertices()
{
}

void
Detector::run()
{
    ros::spin();
}

void
Detector::camera_cb(const image_p& depth_msg,
		    const camera_info_p& camera_info_msg)
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

void
Detector::cloud_cb(const cloud_p& cloud_msg)
{
    try
    {
	using namespace	sensor_msgs;
	using namespace	aist_utility;
	
	_cloud.resize(cloud_msg->height, cloud_msg->width);
	pointcloud_to_points<float>(*cloud_msg, _cloud.begin(),
				    milimeters<float>);
	
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
