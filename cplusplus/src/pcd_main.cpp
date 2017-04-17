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
#include <map>
#include "../include/viz.h"
#include "../include/helper.h"
#include "../include/sampling.h"
#include "../include/segmentation.h"


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

  int bool_viz_or = 0;
  int bool_viz_rs = 0;
  
  if(bool_viz_or){
    printf("visualization of original point cloud! \n");
    pc_viz(pc);
  }


  /************* STAGE 1 : SAMPLING ************/

  /* STEP 1:
     For all i, points:
     In single step, compute (iterative sum/reduction) for all j neighbors of i, A_i,j * x_j, and sum. then || x_i - sum || is weight
     A_i,j is a kernel: RBF/Gaussian/Epachanikov. Store in importance vector. Normalize.
     Simultaneously compute density. (can do this later since only need this for sampled points)
  */
  // tuning_parameters
  float sigma_sq = 0.00005;
  int numNbrs = 20;

  int num_pts = pc->size();
    // compute importance weights
  std::vector<float> imp_wt(num_pts,0.0);
  computeWeights(pc, num_pts, sigma_sq, numNbrs,imp_wt);
  printf("computed weights\n");

  /* STEP 2:
     Randomly Sample K points according to weights , and construct new resampled pointcloud.
   */

// get sample indices
  int total_samples = 10000;
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
  int num_pts_rs = total_samples;
  float sigma_sq_seg = 0.00005;
  int numNbrs_seg = 50;

  std::vector<float> p_density(num_pts_rs,0.0);
  printf("about to compute density");
  segmentation_computeDensity(pc_rs,num_pts_rs,sigma_sq_seg,numNbrs_seg,p_density);
  std::vector<int> parents(num_pts_rs,0);
  std::vector<float> distances(num_pts_rs,0.0);

 //segmentation parameter:
  float tau = 0.1;
  segmentation_linkNeighbors(pc_rs,num_pts_rs, tau, p_density, parents,distances);
  constructSegments(pc_rs,num_pts_rs,parents,distances);

  //print parents:
  int num_trees = 0;
  for (int i =0 ; i<num_pts_rs; i++){
    //printf("i:%d , parent:%d \n",i,parents[i]);
    if(parents[i]==i){
      num_trees++;
    }
  }
  printf("finished segmentation! \n");
  printf("number of segments: %d \n", num_trees);


/************* visualize segments: **********/


//rootlist:
   std::vector<int> clusterids;
  for(int i = 0; i < num_pts_rs ; i++){
    if (parents[i]==i)
    {
      clusterids.push_back(i);
    }
  }

  printf("created clusterIDs \n ");

  std::map< int,std::vector<int> > parent_map;
//construct segments:
  for (std::vector<int>::iterator it = clusterids.begin(); it != clusterids.end(); ++it)
  {
    parent_map.insert(std::pair<int,std::vector<int> >(*it,std::vector<int>()));
  }

  for (int i = 0; i < num_pts; ++i)
  {
    parent_map[parents[i]].push_back(i);
  } 

  printf("created map between clusterIds and members \n ");

  std::map<int,int> red_map;
  std::map<int,int> green_map;
  std::map<int,int> blue_map;

// iterate over each segment:
  for (std::vector<int>::iterator i = clusterids.begin(); i != clusterids.end(); ++i)
  { 
    int red_sum = 0;
    int green_sum = 0;
    int blue_sum = 0;
    int count = 0;
    // iterate over each point in segment
    for (std::vector<int>::iterator j = parent_map[*i].begin(); j != parent_map[*i].end(); ++j)
    {
      int red = static_cast< int >(pc_rs->points[*j].r);
      int green = static_cast< int >(pc_rs->points[*j].g);
      int blue = static_cast< int >(pc_rs->points[*j].b);

      blue_sum += blue;
      green_sum += green;
      red_sum += red;
      count++;
    }

    red_sum =  red_sum/count;
    green_sum =  green_sum/count;
    blue_sum =blue_sum/count;
    red_map[*i] = red_sum;
    green_map[*i] = green_sum;
    blue_map[*i] = blue_sum;
  } 

printf("computed color composition for each segment (avg for viz) \n ");

// construct new point cloud:

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg (new pcl::PointCloud<pcl::PointXYZRGB>);
  pc_seg->width = num_pts_rs;
  pc_seg->height = 1;
  pc_seg->points.resize(pc_rs->height * pc_rs->width);
  int idx = 0;
  for (std::vector<int>::iterator i = clusterids.begin(); i != clusterids.end(); ++i){

    uint8_t r = (uint8_t) red_map[*i];
    uint8_t g = (uint8_t) green_map[*i];
    uint8_t b = (uint8_t) blue_map[*i];
    for (std::vector<int>::iterator j = parent_map[*i].begin(); j != parent_map[*i].end(); ++j){

      int pt_id = *j;
    pc_seg->points[idx].x = pc_rs->points[pt_id].x;
    pc_seg->points[idx].y = pc_rs->points[pt_id].y;
    pc_seg->points[idx].z = pc_rs->points[pt_id].z;
    pc_seg->points[idx].r = r;
    pc_seg->points[idx].g = g;
    pc_seg->points[idx].b = b;
    idx++;
    }
  }


printf("constructed new point cloud\n");



int bool_viz_seg = 1;
printf("visualization of segmented point cloud!\n");
  if(bool_viz_seg){
    pc_viz(pc_seg);
  }



































  return (0);

}
