#include <iostream>
#include <stdio.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/console/parse.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/ply_io.h>
#include <algorithm>
#include <functional>
#include <math.h>       /* pow */
#include <cstdlib>
#include <string>
#include <map>
#include "../include/viz.h"
#include "../include/helper.h"
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh.h>
#include "../include/obj_features.h"
// Types
typedef pcl::PointXYZRGB PointT;

using namespace pcl;

int
main (int argc, char** argv)
{ 


  printf("reading point cloud file! \n");
  //  std::string filename("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/rgbd-scenes_aligned/kitchen_small/kitchen_small_1/kitchen_small_1.pcd");
//  std::string filename("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/sample.pcd");

 std::string filename1("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/objects_3dwarehouse/pc/cereal_box/cereal_box_3_3.ply");
 std::string filename2("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/objects_3dwarehouse/pc/cereal_box/cereal_box_5_0.ply");
 std::string filename3("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/objects_3dwarehouse/pc/soda_can/soda_can_1_0.ply");

 pcl::PointCloud<pcl::PointXYZ>::Ptr pc1 (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::io::loadPLYFile(filename1,*pc1);
 pcl::PointCloud<pcl::PointXYZ>::Ptr pc2 (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::io::loadPLYFile(filename2,*pc2);
   pcl::PointCloud<pcl::PointXYZ>::Ptr pc3 (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::io::loadPLYFile(filename3,*pc3);
printf("done reading! \n");


  int bool_viz_or = 1;  
  if(bool_viz_or){
    printf("visualization of original point cloud! \n");
    pc_viz(pc1);
    pc_viz(pc2);
    pc_viz(pc3);
   }


  //get average histogram:
std::vector<float> avg_hist1(33,0.0);  
computeFeatureHistogram(pc1,avg_hist1);
std::vector<float> avg_hist2(33,0.0);  
computeFeatureHistogram(pc2,avg_hist2);
std::vector<float> avg_hist3(33,0.0);  
computeFeatureHistogram(pc3,avg_hist3);


printf("distance between 1 and 2: %f", computeFeatureHistogramDist(avg_hist1,avg_hist2));
printf("distance between 1 and 3: %f", computeFeatureHistogramDist(avg_hist1,avg_hist3));
















}
