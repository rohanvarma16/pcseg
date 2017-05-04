#include <math.h>      
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector> 
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include "../include/helper.h"
#include <map>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


void constructSegmentedPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg,pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs, int num_pts_rs,std::vector<int>& parents,float* flattenXYZ_rs,float* flattenRGB_rs){

    std::vector<int> clusterids;
    int counter = 0;
    for(int i = 0; i < num_pts_rs ; i++){
      if (parents[i]==i)
      {
	counter++;
        clusterids.push_back(i);
      }
    }

    printf("num_clusters: %d, \n",counter);
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
	int red = (int) floor((flattenRGB_rs[(3 * (*j))] * 255.0));
	int green= (int) floor((flattenRGB_rs[(3 * (*j)) + 1] * 255.0));
	int blue= (int) floor((flattenRGB_rs[(3 * (*j)) + 2] * 255.0));
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

    pc_seg->width = num_pts_rs;
    pc_seg->height = 1;
    pc_seg->points.resize(pc_rs->height * pc_rs->width);
    int idx = 0;
    srand (time(NULL));

    for (std::vector<int>::iterator i = clusterids.begin(); i != clusterids.end(); ++i){

      uint8_t red_comp = (uint8_t) red_map[*i];
      uint8_t green_comp = (uint8_t) green_map[*i];
      uint8_t blue_comp = (uint8_t) blue_map[*i];
      int size_cluster = 0;

      for (std::vector<int>::iterator j = parent_map[*i].begin(); j != parent_map[*i].end(); ++j){

        int pt_id = *j;

	pc_seg->points[idx].x = flattenXYZ_rs[3 * pt_id];
	pc_seg->points[idx].y = flattenXYZ_rs[(3 * pt_id) + 1];
	pc_seg->points[idx].z = flattenXYZ_rs[(3 * pt_id) + 2];
        pc_seg->points[idx].r = red_comp;
        pc_seg->points[idx].g = green_comp;
        pc_seg->points[idx].b = blue_comp;
        idx++;
	size_cluster++;
      }
      //      printf("cluster size: %d, r: %d, g: %d, b: %d \n",size_cluster,red_comp,green_comp,blue_comp);
    }
    printf("constructed segmented point cloud\n");
  }













