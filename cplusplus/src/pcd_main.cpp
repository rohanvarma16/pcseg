#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/console/parse.h>
#include <boost/program_options.hpp>

#include <algorithm>
#include <functional>
#include <math.h>       /* pow */
#include <cstdlib>
#include <string>
#include <map>
//#include "../include/viz.h"
#include "../include/helper.h"
#include "../include/sampling.h"
#include "../include/segmentation.h"
#include "../include/voxel.h"
#include "../include/cycleTimer.h"

// Types
typedef pcl::PointXYZRGB PointT;


void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-n number of points"<<endl;
  cout<<"         "<<"-m number of samples"<<endl;
  cout<<"         "<<"-a voxelization grid size in sampling step"<<endl;
  cout<<"         "<<"-b voxelization grid size in segmentation step"<<endl;
}


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




  int num_samples_in;
  float grid_size_1;
  float grid_size_2;
  int num_pts_in;
  int help_bool = 0;
  int full_file = 1;
  // default:

  if(argc == 1)
    {
      printf("no args!, running default parameters \n");
      num_samples_in = 100000;
      grid_size_1 = 20;
      grid_size_2 = 0.12;
      full_file = 0;
      num_pts_in = 1000000;
      cout<<" number of pts input: "<<num_pts_in<<endl;
      cout<<" number of samples: "<<num_samples_in<<endl;
      cout<<" grid size in sampling: "<<grid_size_1<<endl;
      cout<<" grid size in segmentation: "<<grid_size_2<<endl;
    }
  else{
    char tmp;
    while((tmp=getopt(argc,argv,"hn:m:a:b:"))!=-1)
      {
	switch(tmp)
	  {
	    /*option h : help infomation*/
	  case 'h':
	    showhelpinfo(argv[0]);
	    help_bool = 1;
	    break;
	    /*option n : set number of input points to read */
	  case 'n':
	    num_pts_in = atoi(optarg);
	    full_file = 0;
	    cout<<" number of pts input: "<<optarg<<endl;
	    break;
	    /*option m : set number of samples */
	  case 'm':
	    num_samples_in = atoi(optarg);
	    cout<<" number of samples: "<<optarg<<endl;
	    break;
	    /*option a : set sampling grid size*/
	  case 'a':
	    grid_size_1 = atof(optarg);
	    cout<<" grid size in sampling: "<<optarg<<endl;
	    break;
	    /*option b: set segmentation grid size */
	  case 'b':
	    grid_size_2 = atof(optarg);
	    cout<<" grid size in segmentation: "<<optarg<<endl;
	    break;
	  default:
	    showhelpinfo(argv[0]);
	    break;
	  }
      }
  }
  if(help_bool == 1){
    return(0);
  }


  double readTime = 0.f;
  double samplingTime = 0.f;
  double segmentTime = 0.f;
  double totalTime = 0.f;

  double startTime = CycleTimer::currentSeconds();


  /************************** STAGE 0: INITIALIZATION ****************************************/

  printf("reading point cloud file! \n");
  //  std::string filename("/Users/rohan/Dropbox/Work/workspace/pointcloud/data/rgbd-scenes_aligned/kitchen_small/kitchen_small_1/kitchen_small_1.pcd");
  std::string filename("/afs/andrew.cmu.edu/usr18/rohanv/data/kitchen_small_1.pcd");
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc (new pcl::PointCloud<pcl::PointXYZRGB>);

  if(loadPC(pc,filename)){
    std::cout << "Loaded "
    << pc->width * pc->height
    << " data points from kitchen_small_1.pcd"
    << std::endl;
  }


  if(full_file == 1){
    num_pts_in = pc->size();
    cout<<" number of pts_input: "<<num_pts_in<<endl;
  }
  int bool_viz_or = 0;
  int bool_viz_rs = 0;
  
  if(bool_viz_or){
    printf("visualization of original point cloud! \n");
    //  pc_viz(pc);
  }


  printf("finished reading point cloud file! \n");
  double endReadTime = CycleTimer::currentSeconds();

  
  
  /************* STAGE 1 : SAMPLING ************/

  /* STEP 1:
     For all i, points:
     In single step, compute (iterative sum/reduction) for all j neighbors of i, A_i,j * x_j, and sum. then || x_i - sum || is weight
     A_i,j is a kernel: RBF/Gaussian/Epachanikov. Store in importance vector. Normalize.
     Simultaneously compute density. (can do this later since only need this for sampled points)
  */
  // tuning_parameters
  

  float sigma_sq = 0.00005;
  //  int numNbrs = 50;
  int numNbrs = (int) grid_size_1;
  int num_pts = num_pts_in;
    // compute importance weights
  std::vector<float> imp_wt(num_pts,0.0);
  computeWeights(pc, num_pts, sigma_sq, numNbrs,imp_wt);
  printf("computed weights\n");

  /* STEP 2:
     Randomly Sample K points according to weights , and construct new resampled pointcloud.
   */

// get sample indices

  int total_samples = num_samples_in;
  std::vector<int> sampleIdx = weightedRandomSample(imp_wt,num_pts,total_samples);
  printf("randomly sampled %d points \n" ,total_samples);

// construct resampled point cloud from sample indices:
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs (new pcl::PointCloud<pcl::PointXYZRGB>);
  resamplePC(pc,pc_rs,sampleIdx,total_samples);
  printf("constructed resampled point cloud \n");

  printf("visualization of resampled point cloud!\n");
  if(bool_viz_rs){
    //    pc_viz(pc_rs);
  }

  double endSamplingTime = CycleTimer::currentSeconds();

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
  segmentation_computeDensity(pc_rs,num_pts_rs,sigma_sq_seg,numNbrs_seg,p_density);
    printf("computed density \n");

  std::vector<int> parents(num_pts_rs,0);
  std::vector<float> distances(num_pts_rs,0.0);

 //segmentation parameter:
float  tau = grid_size_2;
  //  float tau = 0.1;
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

  double endTime = CycleTimer::currentSeconds();

  readTime               = 1000.f * (endReadTime - startTime);
  samplingTime           = 1000.f * (endSamplingTime - endReadTime);
  segmentTime            = 1000.f * (endTime - endSamplingTime);
  totalTime              = 1000.f * (endTime - startTime);

  printf("Time for reading input PCD file:              %.4f ms\n", readTime);
  printf("Time for sampling:                            %.4f ms\n", samplingTime);
  printf("Time for segmentation:                        %.4f ms\n", segmentTime);
  printf("Total time:                                   %.4f ms\n", totalTime);



/************* visualize segments: **********/

// construct segmented point cloud:

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg (new pcl::PointCloud<pcl::PointXYZRGB>);
  constructSegmentedPC(pc_seg,pc_rs,num_pts_rs,parents);

  int bool_viz_seg = 0;
  printf("visualization of segmented point cloud!\n");
  if(bool_viz_seg){
    //  pc_viz(pc_seg);
  }

  return (0);

}
