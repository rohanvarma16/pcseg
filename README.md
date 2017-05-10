# Parallel Point Cloud Processing and Segmentation
## Ardra Singh (ardras)
## Rohan Varma (rohanv)

<div style="text-align: center;"><a class="nav"  href="https://rohanvarma16.github.io/pcseg/proposal" target="_blank">Proposal</a></div>

<div style="text-align: center;"><a class="nav"  href="https://rohanvarma16.github.io/pcseg/checkpoint" target="_blank">Checkpoint Report</a></div>

### Introduction and Brief Overview:

In our project, we study spatially local computation on 3-d pointclouds of sizes in the millions and present two main ways with which we can increase speed/throughput on a GPU while preserving the performance accuracy of the end result. By spatially-local, we mean that each pixel independently, performs a computation based on the points in its local neighborhood.
There are two main contributions we want to highlight here. First, we show a mapping of approximate local neighborhoods to CUDA thread blocks to accelerate GPU throughput while preserving accuracy.
Secondly, we implement a fast parallel version of the contour-preserving resampler presented in [this paper] to subsample the pointcloud (preserving just 5% of the points does well!) while preserving important features. This helps to alleviate the cost of the high redundancy, while still preserving output accuracy.

Specifically we chose to study the important problem of segmentation which is an important step in many computer vision application pipelines. That is,  “clustering” point clouds into multiple homogeneous regions, the points in the same region will have the same properties. This problem has many applications in robotics such as intelligent vehicles, autonomous mapping, navigation, household tasks and so on. Point clouds are a particular challenge because they often have uneven point density, high redundancy, and are noisy with many outliers. 

Since it not clear what is the best “metric” for evaluating segmentation, we build a simple object detector (by computing features and finding the nearest neighbor in a pre-computed object feature database). On a high level, we say a segmentation of good quality, if it’s able to reliably detect objects which we have trained upon. We would like to emphasize here, that the object detection is not the main focus of our project. Instead, we show how we can perform fast processing of point-cloud data on a GPU.

We primarily use the <div style="text-align: center;"><a class="nav"  href="http://rgbd-dataset.cs.washington.edu/dataset/" target="_blank">RGB-D Object Database</a></div> which has point clouds of scenes with objects like below as well as models of the objects themselves which we use to train features (using the point cloud library).
<img src="pc_or.png">

### Design and Challenges:

## v 1.0
<img src="block1.png">

The block diagram of the initial design is as above. We use the quick shift algorithm to perform the image segmentation. It is amenable to parallelism because of it’s computational characteristics and memory access patterns.
As can be seen in the block diagram, we first compute a local density estimate of each point by looking at it’s neighborhood, before again doing a spatially local computation to construct a tree before cutting the tree appropriate to get the resulting segmentation.

The first contribution is to leverage the spatial local characteristics of the computation to voxelize the pointcloud and map each voxel to a CUDA thread block. This way every point in a voxel performs the same computation over its neighborhood and possesses the same memory access patterns. This change to the original framework makes it well suited for a fast CUDA implementation. We voxelize by cubing the minimum bounding box of the point cloud. The neighborhood of a voxel is it's neighboring voxels.
<img src="voxelgrid.png">
<img src="voxel_nbr.jpg">


However, this is quite slow, especially because the large number of points being processed in the segmentation and tree-cutting step. We ask whether we can do better? 


## v2.0

<img src="block2.png">
Here we introduce our second main contribution, a resampling block into our design. The sampling step occurs in two stages, we first need to assign an importance weight to each point (effectively a local high pass filter, again a spatially local computation), before performing a weighted sampling of the points. The latter can be implemented quickly with the help of CUDA  thrust libraries. It turns out we can subsample by preserving up to 5% of the total number of points and preserve detection performance.

We mention a caveat here in that, our weighted CPU-based sequential sampler, performs the sampling naively with O(KN) (K is the number of samples). (The sampler pre-computes a rolling sum of the normalized weights and then samples a uniform random number and see in which bin it falls (binary search)). While we implement this same algorithm on CUDA, this is not the fastest way to perform weighted sampling on CUDA. A faster CPU-based implementation would be based on the Alias-Walker method which samples in O(K+N). 
Another interesting point is that the sampler helps to smoothen the uneven density of points across the space. 

### Results:
Below is an image of the result of the segmentation on the kitchen scene. The original point cloud has around 3 million points and we preserve only 80000 samples.
<img src="seg_example.png">
Example of a a cereal box we detect as a result of the segmentation.
<img src="detector.png">

Below, we present preliminary analysis of our results.
## GPU vs. CPU:
We provide a caveat here that this is with respect to a single-threaded sequential implementation. For our final results, we intend to compare the GPU version with a multi-threaded version.

## Time with Sampling vs. Without Sampling:

We note that the overhead of sampling is largely negligible compared to the drastic speedup of the segmentation block. We note that we sample ?? points which preserves detection performance. This is for the CUDA-based implementation.

## Comparison with Segmentation of an Image:

Images exhibit regularity and here we intend to analyze how much the irregularity of the density points in the space affects our performance. Irregularity affects the number of points each block/voxel processes. This in turn leads to extremely skewed workload imbalance patterns. We equalize the number of pixels in the image and the number of points in the point cloud. 

### Summary
This is simply a high level overview of our project and design. We have tried to highlight the novel and most interesting parts of our design. For our final report/presentation, we intend to show more detailed performance analysis of each step in the design flow and a more comprehensive analysis of the workload, the introduction of the sampling stage, and work on further tuning the object detector.











