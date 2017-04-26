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

#include "../include/helper.h"
#include "../include/sampling.h"
#include "../include/segmentation.h"
#include <omp.h>
#include "../include/cuda_test.h"


 

// Types
typedef pcl::PointXYZRGB PointT;


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
  printf("reading point cloud file! \n");
  //std::string filename("/afs/andrew.cmu.edu/usr18/ardras/private/15-618/pcseg/cuda/data/sample.pcd");
  std::string filename("/afs/andrew.cmu.edu/usr18/ardras/data/kitchen_small_1.pcd");

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc (new pcl::PointCloud<pcl::PointXYZRGB>);

  if(loadPC(pc,filename)){
    std::cout << "Loaded "
	      << pc->width * pc->height
	      << " data points from kitchen_small_1.pcd"
	      << std::endl;
  }



  int num_pts = pc->size();
  //num_pts = 1000000;
  /**** pre-processing *****/
  float min_x = pc->points[0].x;
  float min_y =pc->points[0].y;
  float min_z = pc->points[0].z; 
  float max_x =pc->points[0].x;
  float max_y =pc->points[0].y;
  float max_z = pc->points[0].z;

  for (int i =1 ; i < num_pts;i++){
    if(pc->points[i].x< min_x)
      min_x = pc->points[i].x;
    if(pc->points[i].y< min_y)
      min_y = pc->points[i].y;
    if(pc->points[i].z< min_z)
      min_z = pc->points[i].z;
    if(pc->points[i].x > max_x)
      max_x = pc->points[i].x;
    if(pc->points[i].y > max_y)
      max_y = pc->points[i].y;
    if(pc->points[i].z > max_z)
      max_z = pc->points[i].z;
  }

  printf("min_x  %0.4f\n", min_x);
  printf("min_y %0.4f\n ",min_y);
  printf("min_z %0.4f \n",min_z);
  printf("max_x  %0.4f\n", max_x);
  printf("max_y %0.4f\n ",max_y);
  printf("max_z %0.4f \n",max_z);

  // initialize voxels

  float x_grid = 0.1;
  float y_grid = 0.1;
  float z_grid = 0.1;
  float x_range = max_x - min_x;
  float y_range = max_y - min_y;
  float z_range = max_z - min_z;
    
  //float max_xyz = std::max(max_x-min_x,max_y-min_y);
  //max_xyz = std::max(max_xyz, max_z-min_z);
  //printf(" = %f\n", max_xyz);
  
  printf("x_range = %f, y_range = %f, z_range = %f\n", x_range, y_range, z_range);

  int x_idx = std::max(1.0, ceil(x_range/x_grid));
  int y_idx  = std::max(1.0, ceil(y_range/y_grid));
  int z_idx = std::max(1.0, ceil(z_range/z_grid));
  int yz_idx = y_idx * z_idx;

  printf("xidx: %d \n" ,x_idx);
  printf("yidx: %d \n" ,y_idx);
  printf("zidx: %d \n" ,z_idx);
  
  int num_voxels = x_idx * y_idx * z_idx;
  printf("num_voxels %d \n", num_voxels);

  std::vector<std::vector<float> > voxel_pointXYZ;
  std::vector<std::vector<float> > voxel_pointRGB;

  for(int i = 0 ; i < num_voxels ; i++){
    voxel_pointXYZ.push_back(std::vector<float>());
    voxel_pointRGB.push_back(std::vector<float>());
  }
    
  for(int i = 0 ; i < num_pts ; i++){
    int voxel_x = (int) floor((pc->points[i].x - min_x)/x_grid);
    int voxel_y = (int) floor((pc->points[i].y - min_y)/y_grid);
    int voxel_z = (int) floor((pc->points[i].z - min_z)/z_grid);    
    int voxel_id = yz_idx * voxel_x + z_idx * voxel_y + voxel_z ;
    
    if(voxel_id > num_voxels){
      std::cout << "DANGER" << endl;
      printf("danger: %0.4f,%0.4f,%0.4f\n",pc->points[i].x,pc->points[i].y,pc->points[i].z);
      printf("danger: %d, %d, %d, %d\n",voxel_x,voxel_y,voxel_z,voxel_id);
    }
    else{
      //    printf("xyz: %0.4f,%0.4f,%0.4f\n",pc->points[i].x,pc->points[i].y,pc->points[i].z );
    voxel_pointXYZ[voxel_id].push_back(pc->points[i].x);
    voxel_pointXYZ[voxel_id].push_back(pc->points[i].y);
    voxel_pointXYZ[voxel_id].push_back(pc->points[i].z);
    voxel_pointRGB[voxel_id].push_back(float(pc->points[i].r)/255.0);
    voxel_pointRGB[voxel_id].push_back(float(pc->points[i].g)/255.0);
    voxel_pointRGB[voxel_id].push_back(float(pc->points[i].b)/255.0);
    }
  }
  
  /*printf("Voxel 32 first point x = %f, y %f, z %f\n", voxel_pointXYZ[32][0], 
                                 voxel_pointXYZ[32][1], voxel_pointXYZ[32][2]);
  printf("Voxel 32 first point r = %f, g %f, b %f\n", voxel_pointRGB[32][0], 
                                 voxel_pointRGB[32][1], voxel_pointRGB[32][2]);
  printf("total size:%d\n", voxel_pointXYZ.size());
  */

  float* flattenXYZ = (float*) malloc(num_pts * 3 * sizeof(float));
  float* flattenRGB = (float*) malloc(num_pts * 3 * sizeof(float));
  int offset_rs = 0;
  int* voxel_offset = (int*) malloc( sizeof(int) * (num_voxels+1));

  for (int i = 0 ; i < num_voxels ; i++){
    voxel_offset[i] = offset_rs;
    offset_rs += (voxel_pointXYZ[i].size() / 3);
    //printf("voxel %d has %d points\n", i, (voxel_pointXYZ[i].size()/3));
  }
  voxel_offset[num_voxels] = offset_rs;
  printf("final offset: %d \n" , offset_rs);
  
  for(int i = 0 ; i < num_voxels ; i++){
    std::copy(voxel_pointXYZ[i].begin(),voxel_pointXYZ[i].end(),flattenXYZ+(3*voxel_offset[i]));
    std::copy(voxel_pointRGB[i].begin(),voxel_pointRGB[i].end(),flattenRGB+(3*voxel_offset[i]));
  }
  
  int* neighbor_ids = (int*) malloc(sizeof(int) * 7 * num_voxels);
  printf("x_idx: %d, y_idx: %d, z_idx:%d, yz_idx: %d \n", x_idx, y_idx, z_idx, yz_idx);
  for(int i = 0 ; i < x_idx ; i++){
    for(int j=0; j < y_idx ; j++){
      for(int k = 0 ; k < z_idx ; k++){
        int voxel_id = yz_idx * i + z_idx * j + k;
	    if(voxel_id >= num_voxels){
	      printf("voxel_id : %d \n", voxel_id);
	      printf("DANGER\n");
	    }
	
	    neighbor_ids[7 * voxel_id] = voxel_id;
	
	    int nbr_id_1 = yz_idx * (i+1) + z_idx * j + k;
	    //if(nbr_id_1 < 0 || nbr_id_1 >= num_voxels)
	    if(i+1 >= x_idx)
          nbr_id_1 = -1;
	    int nbr_id_2 = yz_idx *(i-1) + z_idx * j + k;
	    //if(nbr_id_2 < 0|| nbr_id_2 >= num_voxels)
        if(i-1 < 0)
          nbr_id_2 = -1;
	    int nbr_id_3 = yz_idx * i + z_idx * (j+1) + k;
	    //if(nbr_id_3 < 0|| nbr_id_3 >= num_voxels)
        if(j+1 >= y_idx)
          nbr_id_3 = -1;
	    int nbr_id_4 = yz_idx * i + z_idx * (j-1) + k;
	    //if(nbr_id_4 < 0|| nbr_id_4 >= num_voxels)
        if(j-1 < 0)
          nbr_id_4 = -1;
	    int nbr_id_5 = yz_idx * i + z_idx * j + (k+1);
	    //if(nbr_id_5 < 0|| nbr_id_5 >= num_voxels)
        if(k+1 >= z_idx)  
          nbr_id_5 = -1;
	    int nbr_id_6 = yz_idx *i + z_idx * j + (k-1);
	    //if(nbr_id_6 < 0|| nbr_id_6 >= num_voxels)
        if(k-1 < 0)  
          nbr_id_6 = -1;
	
	    neighbor_ids[7*voxel_id + 1] = nbr_id_1;
	    neighbor_ids[7*voxel_id+ 2] = nbr_id_2;
	    neighbor_ids[7*voxel_id+ 3] = nbr_id_3;
	    neighbor_ids[7*voxel_id+ 4] = nbr_id_4;
	    neighbor_ids[7*voxel_id+ 5] = nbr_id_5;
	    neighbor_ids[7*voxel_id+ 6] = nbr_id_6;
        //printf("voxel %d has neighbours %d %d %d %d %d %d\n", voxel_id, nbr_id_1,
        //        nbr_id_2, nbr_id_3, nbr_id_4, nbr_id_5, nbr_id_6);
      }
    }
  }

  /*
  for(int i = 0 ;i < num_voxels ; i++){
    printf("voxel_id: %d, mynbrids: ",i);
    for(int j = 0; j < 7 ; j++){
      printf(" %d , ",neighbor_ids[7*i + j]);
    }
    printf("\n");

  }
  */

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


    
    int num_samples = 20000;
    uint* samples_arr = (uint*)malloc(num_samples* sizeof(uint));
    device_setup(num_pts, num_voxels, flattenXYZ,flattenRGB,voxel_offset,
               neighbor_ids,x_idx,y_idx,z_idx,num_samples,samples_arr);


    for(int i = 0 ; i < 10; i++){
        printf("samples_arr[%d] = %d , " , i ,samples_arr[i] );
    }
  

    // construct resampled point cloud from sample indices:
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_rs (new pcl::PointCloud<pcl::PointXYZRGB>);
    //    resamplePC(pc,pc_rs,samples_arr,num_samples);
  
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
    std::cerr << "Saved " << pc_rs->size () << " data points to rs_pcd.pcd." << std::endl;


    printf("constructed and wrote resampled point cloud \n");
    //

    //           STAGE 2 : SEGMENTATION              //
    // STEP 3:
    //   Compute P_n = sum A'_{i,j} for all i in resampled point cloud. A'_{i,j} is another kernel.

    // STEP 4:
    // Do segmentation procedure for each i and link to highest parent. 
    //
    int num_pts_rs = num_samples;
    float sigma_sq_seg = 0.00005;
    int numNbrs_seg = 50;

    //segmentation parameter:
    float tau = 0.1;
  
    //print parents:
  
    printf("finished segmentation! \n");
    //  printf("number of segments: %d \n", num_trees);
        

    /************* visualize segments: **********/

    // construct segmented point cloud:

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pc_seg (new pcl::PointCloud<pcl::PointXYZRGB>);
    //constructSegmentedPC(pc_seg,pc_rs,num_pts_rs,parents);
    return (0);

}
