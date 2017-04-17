#include <math.h>      
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector> 
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/kdtree/kdtree_flann.h>
#include "../include/helper.h"
#include <map>

float computeColorDistance(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc,int i,int x_j){
  float d = 0.0;
  d+= pow((float(pc->points[x_j].r)/255.0)-(float(pc->points[i].r)/255.0),2.0);
  d+= pow((float(pc->points[x_j].g)/255.0)-(float(pc->points[i].g)/255.0),2.0);
  d+= pow((float(pc->points[x_j].b)/255.0)-(float(pc->points[i].b)/255.0),2.0);
  return(d);
}



void segmentation_computeDensity(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, int num_pts, float sigma_sq, int K, std::vector<float>& p_density){

  pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
  kdtree.setInputCloud (pc);
  printf("loaded point cloud into kdtree object for density computation \n");

for(int i = 0 ; i < num_pts ; i++){

    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    pcl::PointXYZRGB searchPoint = pc->points[i];
    float Aij_ew_p;
    float sum_p = 0.0;

    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
      for (size_t j = 1; j < pointIdxNKNSearch.size(); ++j)
      {        
       int x_j = pointIdxNKNSearch[j];
       float pointColorDistance = computeColorDistance(pc,i,x_j);
       Aij_ew_p = exp(-1.0*(pointNKNSquaredDistance[j] +pointColorDistance )/sigma_sq);
       sum_p += Aij_ew_p;
     }
     p_density[i] = sum_p;
   }
   else{
    p_density[i] = 0.0;
  }
}
}

void segmentation_linkNeighbors(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc,int num_pts, float radius,std::vector<float>& p_density, std::vector<int>& parents, std::vector<float>& distances){

	pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
  kdtree.setInputCloud (pc);
  printf("loaded point cloud into kdtree object to link neighbors \n");
  int x_j;
  for(int i = 0 ; i < num_pts ; i++){

    std::vector<int> pointIdxRadiusSearch;
    std::vector<float> pointRadiusSquaredDistance;
    pcl::PointXYZRGB searchPoint = pc->points[i];
    parents[i] = i;
    // do radius search! 
    if ( kdtree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0 )
    { 
      float min_distance = 10000000;
      
      for (size_t j = 1; j < pointIdxRadiusSearch.size(); ++j)
      {        
       x_j = pointIdxRadiusSearch[j];
       float dist = sqrt(pointRadiusSquaredDistance[j]);
       if(p_density[x_j] > p_density[i] && dist<min_distance){
       	distances[i] = dist;
       	parents[i] = x_j;
       	min_distance = dist;
       }
     }
   }
   else{
   	// has no neighbors: lone tree!
    parents[i] = i;
  }
}
}


void constructSegments(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, int num_pts, std::vector<int>& parents, std::vector<float>&distances){

  std::vector<int> node_done(num_pts,0);
  int delta;
  int done = 0;
  int treedepth= 0;

 
  while(! done)
{ delta =0;

  for(int i =0 ; i< num_pts ; i++){

    if(parents[parents[i]] != parents[i]){
      parents[i] = parents[parents[i]];
      delta++;
    }
  }
  if(delta == 0){
    done = 1;
  }
  treedepth++;
}
printf("treedepth: %d \n",treedepth);


}















