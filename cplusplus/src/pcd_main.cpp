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
#include <math.h>       /* pow */
#include <cstdlib>
#include <string>
#include "../include/viz.h"
#include "../include/helper.h"
#include "../include/sampling.h"


int loadPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, std::string filename){

  if (pcl::io::loadPCDFile<pcl::PointXYZRGB> (filename, *pc) == -1) //* load the file
    {
      PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
      return (-1);
    }
  return (1);
}

int
main (int argc, char** argv)
{ 
  /************************** STAGE 0: INITIALIZATION ****************************************/

  printf("reading point cloud file! \n");
  std::string filename("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/rgbd-scenes_aligned/kitchen_small/kitchen_small_1/kitchen_small_1.pcd");
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc (new pcl::PointCloud<pcl::PointXYZRGB>);

  if(loadPC(pc,filename)){
      std::cout << "Loaded "
            << pc->width * pc->height
            << " data points from kitchen_small_1.pcd"
            << std::endl;
  }

  int bool_viz_or = 1;
  int bool_viz_rs = 1;
  
  if(bool_viz_or){
    printf("visualization of original point cloud! \n");
    pc_viz(pc);
  }



  /************* STAGE 1 : SAMPLING ************/

  /* STEP 1:
     For all i, points:
     In single step, compute (iterative sum/reduction) for all j neighbors of i, A_i,j * x_j, and sum. then || x_i - sum ||
     A_i,j is a kernel: RBF/Gaussian/Epachanikov. Store in importance vector. Normalize.
  */
  // tuning_parameters
  float sigma_sq = 0.00005;
  int numNbrs = 20;

  int num_pts = pc->size();
    // compute importance weights
  std::vector<float> imp_wt = computeWeights(pc, num_pts, sigma_sq, numNbrs);
  printf("computed weights\n");

  /* STEP 2:
     Randomly Sample K points according to weights , and construct new resampled pointcloud.
   */

// get sample indices
int total_samples = 1000;
std::vector<int> sampleIdx = weightedRandomSample(imp_wt,num_pts,total_samples);
printf("randomly sampled %d points \n" ,total_samples);

// construct resampled point cloud from sample indices:
pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs (new pcl::PointCloud<pcl::PointXYZRGB>);
resamplePC(pc,pc_rs,sampleIdx,total_samples);
printf("constructed resampled point cloud \n");

printf("visualization of resampled point cloud!\n");
  if(bool_viz_rs){
    pc_viz(pc_rs);
  }

  /*********** STAGE 2 : SEGMENTATION *************/
   // STEP 3:
   //   Compute P_n = sum A'_{i,j} for all i in resampled point cloud. A'_{i,j} is another kernel.
   
  /* STEP 4:
     Do segmentation procedure for each i and link to highest parent. 
   */
//float* distance =
  return (0);

}
