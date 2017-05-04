#ifndef PREPROCESS_H
#define PREPROCESS_H

void preprocess_step1(int num_pts, pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, float* pdensity, float x_grid,float y_grid,
                      float z_grid, int& num_voxels,int& x_idx, int& y_idx, int& z_idx, int& yz_idx,
                      std::vector<std::vector<float> >& voxel_pointXYZ, std::vector<std::vector<float> >& voxel_pointRGB,
                      std::vector<std::vector<float> >& voxel_pdensity);

void preprocess_step2(std::vector<std::vector<float> >& voxel_pointXYZ,std::vector<std::vector<float> >& voxel_pointRGB,
                      std::vector<std::vector<float> >& voxel_pdensity, int x_idx,int y_idx,int z_idx,int yz_idx,int num_voxels,
                      float* flattenXYZ,float* flattenRGB,float* pdensity_new,int* voxel_offset,int* neighbor_ids);
#endif
