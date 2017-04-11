#ifndef VIZ_H
#define VIZ_H

boost::shared_ptr<pcl::visualization::PCLVisualizer> simpleVis (pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud);
boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis (pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr cloud);
void pc_viz (pcl::PointCloud<pcl::PointXYZRGB>::Ptr ptCloud);

#endif