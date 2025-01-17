/*!
* \file		Detector.h
* \author	Toshio UESHIBA
*/
#include <image_transport/image_transport.h>
#include <image_transport/subscriber_filter.h>
#include <sensor_msgs/image_encodings.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>
#include <ddynamic_reconfigure/ddynamic_reconfigure.h>
#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <sensor_msgs/image_encodings.h>
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
    template <class T>
    using point3_t	= cv::Vec<T, 3>;

    template <class T>
    class PointCloud
    {
      public:
	void	resize(int h, int w)
		{
		    _h = h;
		    _w = w;
		    _vertices.resize(_h * _w);
		}
	auto	height()		const	{ return _h; }
	auto	width()			const	{ return _w; }
	bool	get(int v, int u, double& x, double& y, double& z) const
		{
		    const auto	i = v * _w + u;
		    z = _vertices[i][2];
		    if (z == 0 || std::isnan(z))
			return false;
		    x = _vertices[i][0];
		    y = _vertices[i][1];
		    return true;
		}
	auto	begin()				{ return _vertices.begin(); }

      private:
	int				_h, _w;
	std::vector<point3_t<T> >	_vertices; // 3D vertices
    };

    using cloud_t	= PointCloud<float>;

  public:
		Detector(const ros::NodeHandle& nh)			;

    void	run()							;

  private:
    void	detect_plane_cb(const camera_info_p& camera_info_msg,
				const image_p&	     image_msg,
				const image_p&	     depth_msg)		;
    template <class T> cv::Vec<T, 3>
		view_vector(T u, T v)				const	;
    template <class T> cv::Vec<T, 3>
		at(const image_t& depth_msg, int u, int v)	const	;
    template <class T> cv::Vec<T, 3>
		at(const image_t& depth_msg, T u, T v)		const	;

  private:
    ros::NodeHandle					_nh;

  // input camera_info/image stuff
    image_transport::ImageTransport			_it;
    message_filters::Subscriber<camera_info_t>		_camera_info_sub;
    image_transport::SubscriberFilter			_image_sub;
    image_transport::SubscriberFilter			_depth_sub;
    message_filters::Synchronizer<sync_policy_t>	_sync;

  // output stuff
    const image_transport::Publisher			_image_pub;
    const ros::Publisher				_pose_pub;

    ddynamic_reconfigure::DDynamicReconfigure		_ddr;

    double						_planarityTolerance;

    cloud_t						_cloud;
    cv_bridge::CvImage					_seg_img;

    ahc::PlaneFitter<cloud_t>				_plane_fitter;
    std::vector<std::vector<int> >			_plane_vertices;
};
}	// namepsace plane_detector
