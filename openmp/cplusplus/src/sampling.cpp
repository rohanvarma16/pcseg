#include <math.h>      
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector> 
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/kdtree/kdtree_flann.h>
#include "../include/helper.h"
#include <omp.h>

int computeWeights(  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, int num_pts, float sigma_sq, int K,std::vector<float>& imp_wt){

  pcl::KdTreeFLANN<pcl::PointXYZRGB> kdtree;
  kdtree.setInputCloud (pc);
  printf("loaded point cloud into kdtree object for fast neighbor search \n");
  
  float feature_length = 6;
  float normalization_sum = 0.0;


  #pragma omp parallel \
    shared(pc,num_pts,sigma_sq,imp_wt,kdtree,feature_length, normalization_sum)
{
#pragma omp for reduction(+:normalization_sum)
  for(int i = 0 ; i < num_pts ; i++){
    
    std::vector<int> pointIdxNKNSearch(K);
    std::vector<float> pointNKNSquaredDistance(K);
    pcl::PointXYZRGB searchPoint = pc->points[i];
    float Aij_ew;
    imp_wt[i] = 0.0;
    std::vector<float> sum_nbr(feature_length,0.0);

    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
    {
      std::vector<float> x_j_feature(feature_length,0.0);
      for (size_t j = 1; j < pointIdxNKNSearch.size(); ++j)
      {        
       int x_j = pointIdxNKNSearch[j];
       
       Aij_ew = exp(-1.0*(pointNKNSquaredDistance[j]/sigma_sq));

       x_j_feature[0] =  Aij_ew * pc->points[x_j].x; 
       x_j_feature[1] =  Aij_ew * pc->points[x_j].y;
       x_j_feature[2] =  Aij_ew * pc->points[x_j].z;
       x_j_feature[3] = Aij_ew * ((float(pc->points[x_j].r))/255.0);
       x_j_feature[4] = Aij_ew * ((float(pc->points[x_j].g))/255.0);
       x_j_feature[5] = Aij_ew * ((float(pc->points[x_j].b))/255.0);

       sum_nbr = operator+(sum_nbr,x_j_feature);

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
    imp_wt[i] = 0.0;
  }
  normalization_sum += imp_wt[i];  
  }

#pragma omp barrier 
#pragma omp for
for(int i = 0 ; i < num_pts ; i++){
  imp_wt[i] = imp_wt[i] / normalization_sum;
}

 }
return(1);
}


std::vector<int> weightedRandomSample (std::vector<float> weights, int num_pts, int total_samples){

	srand (static_cast <unsigned> (time(0)));
	std::vector<int> samples(num_pts,0);
	std::vector<float> imp_wt_rs(num_pts,0.0);
#pragma omp parallel \
 shared(weights,num_pts,total_samples,samples,imp_wt_rs)
    { 
	// calculate rolling sum of weights:
#pragma omp for	
for (int i = 1 ; i< num_pts ; i++){
   imp_wt_rs[i] = imp_wt_rs[i-1] + weights[i-1];
 }

#pragma omp barrier 

#pragma omp for
 for (int i = 0 ; i < total_samples ; i++){
   float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	  // find bin:
   for (int j = 0 ; j < num_pts-1; j++){
    if(r <= imp_wt_rs[j+1] && r >= imp_wt_rs[j]){
     samples[i] = j;
	        //printf("r: %f, low: %f, high: %f \n",r,imp_wt_rs[j],imp_wt_rs[j+1]);
	        //printf("sampleidx: %d, imp_wt: %0.9f \n ", j, imp_wt[j]);
     break;
   }
 }
 if(r>imp_wt_rs[num_pts-1]){
   samples[i] = num_pts-1;
 }
}
    }

return(samples);
}

void resamplePC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs, 
	std::vector<int> samples, int total_samples){
  pc_rs->width = total_samples;
  pc_rs->height = 1;
  pc_rs->points.resize(pc_rs->height * pc_rs->width);
  #pragma omp parallel for 
  for (int i =0 ; i < total_samples ;i++){
    pc_rs->points[i].x = pc->points[samples[i]].x;
    pc_rs->points[i].y = pc->points[samples[i]].y;
    pc_rs->points[i].z = pc->points[samples[i]].z;
    pc_rs->points[i].r = pc->points[samples[i]].r;
    pc_rs->points[i].g = pc->points[samples[i]].g;
    pc_rs->points[i].b = pc->points[samples[i]].b;
  }
}
