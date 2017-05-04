#ifndef SEGMENTATION_H
#define SEGMENTATION_H
void constructSegmentedPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg,pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs,
			  int num_pts_rs,std::vector<int>& parents,float* flattenXYZ_rs,float* flattenRGB_rs);
#endif
