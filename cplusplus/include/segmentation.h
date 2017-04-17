#ifndef SEGMENTATION_H
#define SEGMENTATION_H

void segmentation_computeDensity(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, int num_pts, float sigma_sq, int K, std::vector<float>& p_density);

void segmentation_linkNeighbors(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc,int num_pts, float radius,std::vector<float>& p_density, std::vector<int>& parents, std::vector<float>& distances);

void constructSegments(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, int num_pts, std::vector<int>& parents, std::vector<float>&distances);
#endif