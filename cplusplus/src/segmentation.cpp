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

  void constructSegmentedPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg,pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs, int num_pts_rs,std::vector<int>& parents){

    std::vector<int> clusterids;
    for(int i = 0; i < num_pts_rs ; i++){
      if (parents[i]==i)
      {
        clusterids.push_back(i);
      }
    }

    printf("created clusterIDs \n");

    std::map< int,std::vector<int> > parent_map;
//construct segments:
    for (std::vector<int>::iterator it = clusterids.begin(); it != clusterids.end(); ++it)
    {
      parent_map.insert(std::pair<int,std::vector<int> >(*it,std::vector<int>()));
    }

    for (int i = 0; i < num_pts_rs; ++i)
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
    printf("constructed segmented point cloud\n");
  }













