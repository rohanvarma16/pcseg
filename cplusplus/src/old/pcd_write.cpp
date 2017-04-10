#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>

int  tester (int x){
  return x+1;
}

int
main (int argc, char** argv)
{

  /* WRITE PCD FILE */
  pcl::PointCloud<pcl::PointXYZ> cloud;

  // Fill in the cloud data
  cloud.width    = 5;
  cloud.height   = 1;
  cloud.is_dense = false;
  cloud.points.resize (cloud.width * cloud.height);

  for (size_t i = 0; i < cloud.points.size (); ++i)
    {
      cloud.points[i].x = 1024 * rand () / (RAND_MAX + 1.0f);
      cloud.points[i].y = 1024 * rand () / (RAND_MAX + 1.0f);
      cloud.points[i].z = 1024 * rand () / (RAND_MAX + 1.0f);
    }

  pcl::io::savePCDFileASCII ("test_pcd.pcd", cloud);
  std::cerr << "Saved " << cloud.points.size () << " data points to test_pcd.pcd." << std::endl;

  for (size_t i = 0; i < cloud.points.size (); ++i)
    std::cerr << "    " << cloud.points[i].x << " " << cloud.points[i].y << " " << cloud.points[i].z << std::endl;


  int y = tester(1);
  printf("y = %d \n", y);

  /* READ PCD FILE */
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud2 (new pcl::PointCloud<pcl::PointXYZRGB>);

  if (pcl::io::loadPCDFile<pcl::PointXYZRGB> ("test_pcd.pcd", *cloud2) == -1) //* load the file
    {
      PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
      return (-1);
    }
  std::cout << "Loaded "
            << cloud2->width * cloud2->height
            << " data points from test_pcd.pcd with the following fields: "
            << std::endl;
  for (size_t i = 0; i < cloud2->points.size (); ++i)
    std::cout << "    " << cloud2->points[i].x
              << " "    << cloud2->points[i].y
              << " "    << cloud2->points[i].z << std::endl;
    
  /* VISUALIZATION 
  pcl::visualization::CloudViewer viewer ("Simple Cloud Viewer");
  viewer.showCloud (cloud2);
  while (!viewer.wasStopped ())
    {
    }
*/
  
  return (0);



}
