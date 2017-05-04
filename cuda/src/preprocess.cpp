#include <iostream>
#include <stdio.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/console/parse.h>
#include <algorithm>
#include <functional>
#include <math.h>
#include <cstdlib>
#include <string>
#include <map>

void preprocess_step1(int num_pts, pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, float* pdensity, float x_grid,float y_grid, 
		      float z_grid, int& num_voxels,int& x_idx, int& y_idx, int& z_idx, int& yz_idx, 
		      std::vector<std::vector<float> >& voxel_pointXYZ, std::vector<std::vector<float> >& voxel_pointRGB,
		      std::vector<std::vector<float> >& voxel_pdensity){

  float min_x = pc->points[0].x;
  float min_y =pc->points[0].y;
  float min_z = pc->points[0].z;
  float max_x =pc->points[0].x;
  float max_y =pc->points[0].y;
  float max_z = pc->points[0].z;

  for (int i =1 ; i < num_pts;i++){
    if(pc->points[i].x< min_x)
      min_x = pc->points[i].x;
    if(pc->points[i].y< min_y)
      min_y = pc->points[i].y;
    if(pc->points[i].z< min_z)
      min_z = pc->points[i].z;
    if(pc->points[i].x > max_x)
      max_x = pc->points[i].x;
    if(pc->points[i].y > max_y)
      max_y = pc->points[i].y;
    if(pc->points[i].z > max_z)
      max_z = pc->points[i].z;
  }

  // initialize voxels
  float x_range = max_x - min_x;
  float y_range = max_y - min_y;
  float z_range = max_z - min_z;

  x_idx = std::max(1.0, ceil(x_range/x_grid));
  y_idx  = std::max(1.0, ceil(y_range/y_grid));
  z_idx = std::max(1.0, ceil(z_range/z_grid));
  yz_idx = y_idx * z_idx;

  num_voxels = x_idx * y_idx * z_idx;
  printf("num_voxels %d \n", num_voxels);

  for(int i = 0 ; i < num_voxels ; i++){
    voxel_pointXYZ.push_back(std::vector<float>());
    voxel_pointRGB.push_back(std::vector<float>());
    voxel_pdensity.push_back(std::vector<float>());
  }

  for(int i = 0 ; i < num_pts ; i++){
    int voxel_x = (int) floor((pc->points[i].x - min_x)/x_grid);
    int voxel_y = (int) floor((pc->points[i].y - min_y)/y_grid);
    int voxel_z = (int) floor((pc->points[i].z - min_z)/z_grid);
    int voxel_id = yz_idx * voxel_x + z_idx * voxel_y + voxel_z ;

    if(voxel_id > num_voxels){
      std::cout << "DANGER" << endl;
    }
    else{
      voxel_pointXYZ[voxel_id].push_back(pc->points[i].x);
      voxel_pointXYZ[voxel_id].push_back(pc->points[i].y);
      voxel_pointXYZ[voxel_id].push_back(pc->points[i].z);
      voxel_pointRGB[voxel_id].push_back(float(pc->points[i].r)/255.0);
      voxel_pointRGB[voxel_id].push_back(float(pc->points[i].g)/255.0);
      voxel_pointRGB[voxel_id].push_back(float(pc->points[i].b)/255.0);
      voxel_pdensity[voxel_id].push_back(pdensity[i]);
    }
  }
}


void preprocess_step2(std::vector<std::vector<float> >& voxel_pointXYZ,std::vector<std::vector<float> >& voxel_pointRGB,
                      std::vector<std::vector<float> >& voxel_pdensity, int x_idx,int y_idx,int z_idx,int yz_idx,int num_voxels,
                      float* flattenXYZ,float* flattenRGB,float* pdensity_new,int* voxel_offset,int* neighbor_ids){

  int offset_rs = 0;
  for (int i = 0 ; i < num_voxels ; i++){
    voxel_offset[i] = offset_rs;
    offset_rs += (voxel_pointXYZ[i].size() / 3);
  }

  voxel_offset[num_voxels] = offset_rs;

  for(int i = 0 ; i < num_voxels ; i++){
    std::copy(voxel_pointXYZ[i].begin(),voxel_pointXYZ[i].end(),flattenXYZ+(3*voxel_offset[i]));
    std::copy(voxel_pointRGB[i].begin(),voxel_pointRGB[i].end(),flattenRGB+(3*voxel_offset[i]));
    std::copy(voxel_pdensity[i].begin(),voxel_pdensity[i].end(),pdensity_new + voxel_offset[i]);
  }

  for(int i = 0 ; i < x_idx ; i++){
    for(int j=0; j < y_idx ; j++){
      for(int k = 0 ; k < z_idx ; k++){
        int voxel_id = yz_idx * i + z_idx * j + k;
        if(voxel_id >= num_voxels){
          printf("voxel_id : %d \n", voxel_id);
          printf("DANGER\n");
        }

        neighbor_ids[7 * voxel_id] = voxel_id;
        int nbr_id_1 = yz_idx * (i+1) + z_idx * j + k;

        if(i+1 >= x_idx)
          nbr_id_1 = -1;

        int nbr_id_2 = yz_idx *(i-1) + z_idx * j + k;
        if(i-1 < 0)
          nbr_id_2 = -1;

        int nbr_id_3 = yz_idx * i + z_idx * (j+1) + k;
        if(j+1 >= y_idx)
          nbr_id_3 = -1;

        int nbr_id_4 = yz_idx * i + z_idx * (j-1) + k;
        if(j-1 < 0)
          nbr_id_4 = -1;

        int nbr_id_5 = yz_idx * i + z_idx * j + (k+1);
        if(k+1 >= z_idx)
          nbr_id_5 = -1;

        int nbr_id_6 = yz_idx *i + z_idx * j + (k-1);
        if(k-1 < 0)
          nbr_id_6 = -1;

        neighbor_ids[7*voxel_id + 1] = nbr_id_1;
        neighbor_ids[7*voxel_id+ 2] = nbr_id_2;
        neighbor_ids[7*voxel_id+ 3] = nbr_id_3;
        neighbor_ids[7*voxel_id+ 4] = nbr_id_4;
        neighbor_ids[7*voxel_id+ 5] = nbr_id_5;
        neighbor_ids[7*voxel_id+ 6] = nbr_id_6;

      }
    }
  }
}

