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
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh.h>
#include "../include/obj_features.h"


void computeFeatureHistogram( pcl::PointCloud<pcl::PointXYZ>::Ptr pc, std::vector<float>& avg_hist){


pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
  ne.setInputCloud (pc);

  // Create an empty kdtree representation, and pass it to the normal estimation object.
  // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
  ne.setSearchMethod (tree);
  // Output datasets
  pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
  // Use all neighbors in a sphere of radius 3cm
  ne.setRadiusSearch (0.03);
  // Compute the features
  ne.compute (*normals);

  // Create the FPFH estimation class, and pass the input dataset+normals to it
  pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
  fpfh.setInputCloud (pc);
  fpfh.setInputNormals (normals);
  // alternatively, if cloud is of tpe PointNormal, do fpfh.setInputNormals (cloud);
  fpfh.setSearchMethod (tree);
  // Output 
  pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs (new pcl::PointCloud<pcl::FPFHSignature33> ());
  // Use all neighbors in a sphere of radius 5cm
  // IMPORTANT: the radius used here has to be larger than the radius used to estimate the surface normals!!!
  fpfh.setRadiusSearch (0.05);
  // Compute the features
  fpfh.compute (*fpfhs);

  //get average histogram:
for(int i = 0 ; i < fpfhs->points.size();i++){
    for(int j= 0 ; j < 33 ; j++){
      avg_hist[j]+= fpfhs->points[i].histogram[j];
    }
  }
for(int j=0;j<33;j++){
  avg_hist[j] = avg_hist[j]/fpfhs->points.size();
  //printf("avg_hist[%d]=%f, ",j,avg_hist[j]);
}   
}


float computeFeatureHistogramDist(std::vector<float> &h1,std::vector<float> &h2){
	float dist = 0.f;

	for(int i = 0; i < 33 ; i++)
	{
		dist += pow(h1[i] - h2[i],2.0);
	}
	dist = sqrt(dist);

	return dist;
}













