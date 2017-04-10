#include <iostream>
#include <stdio.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/console/parse.h>
#include "viz.h"
#include <algorithm>
#include <functional>
#include <math.h>       /* pow */
#include <cstdlib>


template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}



int
main (int argc, char** argv)
{
  
  int bool_viz = 0;
  /* READ PCD FILE */
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc (new pcl::PointCloud<pcl::PointXYZRGB>);

  if (pcl::io::loadPCDFile<pcl::PointXYZRGB> ("data/rgbd-scenes_aligned/kitchen_small/kitchen_small_1/kitchen_small_1.pcd", *pc) == -1) //* load the file
    {
      PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
      return (-1);
    }
  std::cout << "Loaded "
            << pc->width * pc->height
            << " data points from test_pcd.pcd with the following fields: "
            << std::endl;

/*
  for (size_t i = 0; i < 5; ++i){
    std::cout << "    " << pc->points[i].x
              << " "    << pc->points[i].y
              << " "    << pc->points[i].z << std::endl;
  std::cout << "    " << (float)pc->points[i].r
	    << " "    << (float)pc->points[i].g
	    << " "    << (float)pc->points[i].b << std::endl;
  }
*/

  int num_pts = pc->size();
  printf("point cloud of size: %d \n", num_pts);

  /************* STAGE 1 : SAMPLING ************/
  /* STEP 1:
     For all i, points:
     In single step, compute (iterative sum/reduction) for all j neighbors of i, A_i,j * x_j, and sum. then || x_i - sum ||
     A_i,j is a kernel: RBF/Gaussian/Epachanikov. Store in importance vector. Normalize.
  */

  pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
  kdtree.setInputCloud (pc);
  printf("loaded point cloud into kdtree object for fast neighbor search \n");

  num_pts = pc->size();
  //float* imp_wt = new float[num_pts];
  std::vector<float> imp_wt(num_pts,0.0);
  float sigma_sq = 0.00005;
  float feature_length = 6;
  float normalization_sum = 0.0;
  float sum;
  for(int i = 0 ; i < num_pts ; i++){
    
    int K = 20;
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    pcl::PointXYZRGB searchPoint = pc->points[i];
    float Aij_ew;
    imp_wt[i] = 0.0;
    std::vector<float> sum_nbr(feature_length,0.0);
    sum = 0.0;
    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
      std::vector<float> x_j_feature(feature_length,0.0);
      for (size_t j = 1; j < pointIdxNKNSearch.size(); ++j)
      {        

        Aij_ew = exp(-1.0*(pointNKNSquaredDistance[j]/sigma_sq));
        //printf("Aij_ew = %f \n", Aij_ew );
        //printf("pointNKNSquaredDistance[j]: %f\n", pointNKNSquaredDistance[j] );
        int x_j = pointIdxNKNSearch[j];
        
        x_j_feature[0] =  Aij_ew * pc->points[x_j].x; 
        x_j_feature[1] =  Aij_ew * pc->points[x_j].y;
        x_j_feature[2] =  Aij_ew * pc->points[x_j].z;
        x_j_feature[3] = Aij_ew * ((float(pc->points[x_j].r))/255.0);
        x_j_feature[4] = Aij_ew * ((float(pc->points[x_j].g))/255.0);
        x_j_feature[5] = Aij_ew * ((float(pc->points[x_j].b))/255.0);

        sum_nbr = operator+(sum_nbr,x_j_feature);
        sum = x_j_feature[0]+ x_j_feature[1]+x_j_feature[2]+x_j_feature[3]+x_j_feature[4]+x_j_feature[5];
      }

      
    //compute norm:
    float norm_sum = 0.0;
    norm_sum += pow(sum_nbr[0]-searchPoint.x,2.0);
    norm_sum += pow(sum_nbr[1] - searchPoint.y,2.0);
    norm_sum += pow(sum_nbr[2] - searchPoint.z,2.0);
    norm_sum += pow(sum_nbr[3] - ((float) searchPoint.r)/255.0,2.0);
    norm_sum += pow(sum_nbr[4]- ((float) searchPoint.g)/255.0,2.0);
    norm_sum += pow(sum_nbr[5] - ((float) searchPoint.b)/255.0,2.0);
    imp_wt[i] = norm_sum; 
    }
    else{
      printf("hello, found no neighbors\n");
    }
    normalization_sum += imp_wt[i];  
}
  
  
  printf("normalization_sum: %f \n", normalization_sum);

  printf("computed weights\n");
/*
  float minwt = 1.0;
  float maxwt = 0.0;
  int counter = 0;
  float sumwt = 0;
  // manual histogram
  for(int i = 0 ; i < num_pts ; i++){
    imp_wt[i] = imp_wt[i] / normalization_sum;

    if(imp_wt[i]<minwt){
      minwt = imp_wt[i];
    }

    if(imp_wt[i] > maxwt){
      maxwt = imp_wt[i];
    }
    if(imp_wt[i] > 0.0000017 ){
      counter++;
      //printf("high importance: %.9f \n", imp_wt[i]);
    }
    sumwt += imp_wt[i];
  }

printf("finished calculating importance weights! \n");
printf("minwt: %0.9f\n", minwt );
printf("maxwt: %0.9f\n", maxwt );
printf("sumwt: %0.9f \n", sumwt );
printf("counter: %d\n", counter );
*/



  /* STEP 2:
     Randomly Sample these points according to weights (1/0), and add to new resampled pointcloud.
   */

for(int i = 0 ; i < num_pts ; i++){
  imp_wt[i] = imp_wt[i] / normalization_sum;
}

srand (static_cast <unsigned> (time(0)));

std::vector<int> sampleIdx(num_pts,0);
int total_samples = 1000;
int num_samples = 0;
std::vector<float> imp_wt_rs(num_pts,0.0);

// rolling sum of weights:
for (int i = 1 ; i< num_pts ; i++){
  imp_wt_rs[i] = imp_wt_rs[i-1] + imp_wt[i-1];
}


for (int i = 0 ; i < total_samples ; i++){
  //printf("iteration: %d \n", i);
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  // find bin:

  for (int j = 0 ; j < num_pts-1; j++){
     if(r <= imp_wt_rs[j+1] && r >= imp_wt_rs[j]){
        num_samples++;
        sampleIdx[i] = j;
        //printf("r: %f, low: %f, high: %f \n",r,imp_wt_rs[j],imp_wt_rs[j+1]);
        //printf("sampleidx: %d, imp_wt: %0.9f \n ", j, imp_wt[j]);
        break;
    }
  }

  if(r>imp_wt_rs[num_pts-1]){
    num_samples++;
    sampleIdx[i] = num_pts-1;
    printf("hello!!!");
  }
}

printf("randomly sampled %d points \n", num_samples );

// construct resampled point cloud from sampleIdx:
pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs (new pcl::PointCloud<pcl::PointXYZRGB>);
pc_rs->width = total_samples;
pc_rs->height = 1;
pc_rs->points.resize(pc_rs->height * pc_rs->width);

for (int i =0 ; i < total_samples ;i++){
  pc_rs->points[i].x = pc->points[sampleIdx[i]].x;
  pc_rs->points[i].y = pc->points[sampleIdx[i]].y;
  pc_rs->points[i].z = pc->points[sampleIdx[i]].z;
  pc_rs->points[i].r = pc->points[sampleIdx[i]].r;
  pc_rs->points[i].g = pc->points[sampleIdx[i]].g;
  pc_rs->points[i].b = pc->points[sampleIdx[i]].b;
}
printf("constructed resampled point cloud");
printf("visualization!\n");
int bool_viz2 = 1;


  if(bool_viz2){
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
  viewer = rgbVis(pc_rs);
  while (!viewer->wasStopped ())
    {
      viewer->spinOnce (100);
      boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
  }

  /*********** STAGE 2 : SEGMENTATION *************/
   // STEP 3:
   //   Compute P_n = sum A'_{i,j} for all i in resampled point cloud. A'_{i,j} is another kernel.
   
  /* STEP 4:
     Do segmentation procedure for each i and link to highest parent. 
   */
//float* distance =

  if(bool_viz){
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
  viewer = rgbVis(pc);
  while (!viewer->wasStopped ())
    {
      viewer->spinOnce (100);
      boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }

  }

  return (0);



}
