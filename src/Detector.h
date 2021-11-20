/*!
* \file		Detector.h
* \author	Toshio UESHIBA
*/
#include <image_transport/image_transport.h>
#include <image_transport/subscriber_filter.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/PointCloud2.h>
#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <cv_bridge/cv_bridge.h>
#include "peac/AHCPlaneFitter.hpp"

namespace plane_detector
{
/************************************************************************
*  class Detector							*
************************************************************************/
class Detector
{
  private:
    using camera_info_t	= sensor_msgs::CameraInfo;
    using camera_info_p	= sensor_msgs::CameraInfoConstPtr;
    using image_t	= sensor_msgs::Image;
    using image_p	= sensor_msgs::ImageConstPtr;
    using sync_policy_t	= message_filters::sync_policies::
			      ApproximateTime<camera_info_t, image_t, image_t>;
    using cloud_t	= sensor_msgs::PointCloud2;
    using cloud_p	= sensor_msgs::PointCloud2ConstPtr;

    template <class T>
    class PointCloud
    {
      public:
	using value_t	= T;
	using vector3_t	= cv::Vec<T, 3>;

      public:
	void	resize(int h, int w)
		{
		    _h = h;
		    _w = w;
		    _points.resize(_h * _w);
		}
	auto	height()		const	{ return _h; }
	auto	width()			const	{ return _w; }
	const vector3_t&
		get(int v, int u)	const	{ return _points[v*_w + u]; }
	auto	begin()				{ return _points.begin(); }

      private:
	int			_h, _w;
	std::vector<vector3_t>	_points; // 3D vertices
    };

    using my_cloud_t	= PointCloud<double>;

  public:
		Detector(const ros::NodeHandle& nh)			;

    void	run()							;

  private:
    void	camera_cb(const image_p& depth_msg,
			  const camera_info_p& camera_info_msg)		;
    void	cloud_cb(const cloud_p& cloud_msg)			;

  private:
    ros::NodeHandle			_nh;

  // input camera_info/image stuff
    image_transport::ImageTransport	_it;
    image_transport::CameraSubscriber	_camera_sub;
    ros::Subscriber			_cloud_sub;

  // output stuff
    const image_transport::Publisher	_image_pub;

    my_cloud_t				_cloud;
    cv_bridge::CvImage			_seg_img;

    ahc::PlaneFitter<my_cloud_t>	_plane_fitter;
    std::vector<std::vector<int> >	_plane_vertices;
};
}	// namepsace plane_detector
