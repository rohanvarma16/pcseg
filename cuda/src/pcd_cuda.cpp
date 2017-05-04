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
#include <algorithm>
#include <functional>
#include <math.h>      
#include <cstdlib>
#include <string>
#include <map>
#include "../include/helper.h"
#include "../include/segmentation.h"
#include "../include/cuda_methods.h"
#include "../include/preprocess.h"
#include <boost/program_options.hpp>

int loadPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc, std::string filename){

  if (pcl::io::loadPCDFile<pcl::PointXYZRGB> (filename, *pc) == -1) //* load the file
  {
    PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
    return (-1);
  }
  return (1);
}

void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl;
  cout<<"         "<<"-n number of points"<<endl;
  cout<<"         "<<"-m number of samples"<<endl;
  cout<<"         "<<"-a voxelization grid size in sampling step"<<endl;
  cout<<"         "<<"-b voxelization grid size in segmentation step"<<endl;
}


int
main (int argc, char** argv)
{ 

  // top-level parameters:
  
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
      grid_size_1 = 0.02;
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


  printf("reading point cloud file! \n");
  //std::string filename("/afs/andrew.cmu.edu/usr18/ardras/private/15-618/pcseg/cuda/data/sample.pcd");
  //  std::string filename("/afs/andrew.cmu.edu/usr18/ardras/data/kitchen_small_1.pcd");
  std::string filename("/afs/andrew.cmu.edu/usr18/rohanv/data/kitchen_small_1.pcd");

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc (new pcl::PointCloud<pcl::PointXYZRGB>);

  if(loadPC(pc,filename)){
    std::cout << "Loaded "
	      << pc->width * pc->height
	      << " data points from kitchen_small_1.pcd"
	      << std::endl;
  }
  /************* STAGE 0 : PREPROCESS POINT CLOUD **************/
  if(full_file == 1){
    num_pts_in = pc->size();
    cout<<" number of pts_input: "<<num_pts_in<<endl;
  }

  int num_pts = num_pts_in;
  float x_grid = grid_size_1;
  float y_grid = grid_size_1;
  float z_grid = grid_size_1;
  int x_idx;
  int y_idx;
  int z_idx;
  int yz_idx;
  int num_voxels;
  std::vector<std::vector<float> > voxel_pointXYZ;
  std::vector<std::vector<float> > voxel_pointRGB;
  std::vector<std::vector<float> > voxel_pdensity_1;
  float* pdensity_dummy = new float[num_pts]();
  float* pdensity_dummy2 = new float[num_pts]();
  preprocess_step1(num_pts,pc,pdensity_dummy,x_grid,y_grid,z_grid,num_voxels,x_idx,y_idx,z_idx,yz_idx,
		   voxel_pointXYZ,voxel_pointRGB,voxel_pdensity_1);

  float* flattenXYZ = (float*) malloc(num_pts * 3 * sizeof(float));
  float* flattenRGB = (float*) malloc(num_pts * 3 * sizeof(float));
  int* voxel_offset = (int*) malloc( sizeof(int) * (num_voxels+1));
  int* neighbor_ids = (int*) malloc(sizeof(int) * 7 * num_voxels);
  preprocess_step2(voxel_pointXYZ,voxel_pointRGB,voxel_pdensity_1,x_idx,y_idx,z_idx,yz_idx,num_voxels,
		   flattenXYZ,flattenRGB,pdensity_dummy2,voxel_offset,neighbor_ids);

  
  /************* STAGE 1 : SAMPLING ************/
  /* STEP 1:
     For all i, points:
     In single step, compute (iterative sum/reduction) for all j neighbors of i, A_i,j * x_j, and sum. then || x_i - sum || is weight
     A_i,j is a kernel: RBF/Gaussian/Epachanikov. Store in importance vector. Normalize.
     Simultaneously compute density. (can do this later since only need this for sampled points)
  */
  /* STEP 2:
     Randomly Sample K points according to weights , and construct new resampled pointcloud.
  */

    int num_samples = num_samples_in;
    uint* samples_arr = (uint*)malloc(num_samples* sizeof(uint));
    float* pdens = (float*) malloc(num_samples * sizeof(float));
    
    cuda_resampling(num_pts, num_voxels, flattenXYZ,flattenRGB,voxel_offset,
	       neighbor_ids,x_idx,y_idx,z_idx,num_samples,samples_arr,pdens);

    // construct resampled point cloud from sample indices:
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs (new pcl::PointCloud<pcl::PointXYZRGB>);  
    pc_rs->width = num_samples;
    pc_rs->height = 1;
    pc_rs->points.resize(pc_rs->height * pc_rs->width);

    for (int i =0 ; i < num_samples ;i++){
      pc_rs->points[i].x = flattenXYZ[samples_arr[i]*3];
      pc_rs->points[i].y = flattenXYZ[(3*samples_arr[i])+1];
      pc_rs->points[i].z = flattenXYZ[(3*samples_arr[i])+2];
      pc_rs->points[i].r = (uint8_t)floor((flattenRGB[samples_arr[i]*3] * 255.0));
      pc_rs->points[i].g = (uint8_t)floor((flattenRGB[(samples_arr[i]*3)+1] *255.0));
      pc_rs->points[i].b = (uint8_t)floor((flattenRGB[(samples_arr[i]*3)+2] *255.0));  
    }

    // write resampled point cloud pc_rs to pcd file 
    pcl::io::savePCDFileASCII ("rs_pcd.pcd", *pc_rs);
    std::cerr << "Saved " << pc_rs->size () << " data points to resampled point cloud rs_pcd.pcd." << std::endl;
    

    /************* STAGE 3 : PREPROCESS POINT CLOUD **************/
    int num_pts_rs = pc_rs->size();
    float x_grid_rs = grid_size_2;
    float y_grid_rs = grid_size_2;
    float z_grid_rs = grid_size_2;
    int x_idx_rs;
    int y_idx_rs;
    int z_idx_rs;
    int yz_idx_rs;
    int num_voxels_rs;
    std::vector<std::vector<float> > voxel_pointXYZ_rs;
    std::vector<std::vector<float> > voxel_pointRGB_rs;
    std::vector<std::vector<float> > voxel_pdensity_rs;

    
    preprocess_step1(num_pts_rs,pc_rs,pdens,x_grid_rs,y_grid_rs,z_grid_rs,num_voxels_rs,
		     x_idx_rs,y_idx_rs,z_idx_rs,yz_idx_rs,voxel_pointXYZ_rs,voxel_pointRGB_rs,voxel_pdensity_rs);
    
    float* flattenXYZ_rs = (float*) malloc(num_pts_rs * 3 * sizeof(float));
    float* flattenRGB_rs = (float*) malloc(num_pts_rs * 3 * sizeof(float));
    float* pdens_new = (float*) malloc(num_pts_rs * sizeof(float));
    int* voxel_offset_rs = (int*) malloc( sizeof(int) * (num_voxels_rs+1));
    int* neighbor_ids_rs = (int*) malloc(sizeof(int) * 7 * num_voxels_rs);
    int* parents = (int*) malloc(num_pts_rs * sizeof(int));
    preprocess_step2(voxel_pointXYZ_rs,voxel_pointRGB_rs,voxel_pdensity_rs,x_idx_rs,y_idx_rs,z_idx_rs,
		     yz_idx_rs,num_voxels_rs,flattenXYZ_rs,flattenRGB_rs,pdens_new,voxel_offset_rs,neighbor_ids_rs);
    
    /*************** STAGE 4 : SEGMENTATION ***************/
    /* STEP 3:
     * Compute P_n = sum A'_{i,j} for all i in resampled point cloud. A'_{i,j} is another kernel.
     * STEP 4:
     * Do segmentation procedure for each i and link to closest parent with a higher density.
     ******************************************************/

    cuda_segmentation(num_pts_rs, num_voxels_rs,pdens_new, flattenXYZ_rs,flattenRGB_rs,voxel_offset_rs,
			neighbor_ids_rs,x_idx_rs,y_idx_rs,z_idx_rs,parents);


    /************* STAGE 5: VISUALIZATION **********/

    // construct segmented point cloud:
    std::vector<int> parents_vec(num_pts_rs);
    std::copy(parents, parents+num_pts_rs,parents_vec.begin()); 
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg (new pcl::PointCloud<pcl::PointXYZRGB>);
    constructSegmentedPC(pc_seg,pc_rs,num_pts_rs,parents_vec,flattenXYZ_rs,flattenRGB_rs);

    // write segmented point cloud pc_rs to pcd file
    pcl::io::savePCDFileASCII ("seg_pcd.pcd", *pc_seg);
    std::cerr << "Saved " << pc_seg->size() << " data points to segmented point cloud seg_pcd.pcd." << std::endl;
    
    return(1);

}
