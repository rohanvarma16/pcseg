## Parallel Point Cloud Processing and Segmentation
Ardra Singh (ardras)
Rohan Varma (rohanv)

You can use the [editor on GitHub](https://github.com/rohanvarma16/pcseg/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Summary
In our project we intend to implement a framework for segmenting point clouds on GPUs. We plan on modifying algorithms based on sequential implementations to design subsampling and segmentation building blocks such that they are amenable to a fast CUDA implementation. We intend on providing a comparison between a CPU-based parallel implementation and an implementation that runs on a GPU. Since segmentation is a critical first step in many applications that process point clouds, a fast implementation that is still accurate and is amenable to processing consecutive frames of point clouds at a high frame rate is valuable.

### Introduction and Background
With the recent development of 3D sensing technologies, 3D point clouds have become an important and practical representation of 3D objects and surrounding environments in many applications, such as virtual reality, mobile mapping, scanning of historical artifacts, 3D printing and digital elevation models. A large number of 3D points on an object's surface measured by a sensing device are called a 3D point cloud. Other than 3D coordinates, a 3D point cloud may also comprise some attributes, such as color, temperature and texture. 3D point cloud segmentation is the process of
classifying point clouds into multiple homogeneous regions, the
points in the same region will have the same properties. The
segmentation is challenging because of high redundancy, uneven
sampling density, and lack explicit structure of point cloud
data. This problem has many applications in robotics such as
intelligent vehicles, autonomous mapping and navigation.
However, segmenting objects in 3D point clouds is not a trivial task. The point
cloud data are usually noisy, sparse, and unorganized. In addition,
the surface shape can be arbitrary with sharp features and
there is no statistical distribution pattern in the data. 




### Framework:
Typically, the number of points in point clouds are on the order of millions. These are however highly redundant and noisy with many outliers. Hence, there is substantial gain from sampling this point cloud into a fewer number of points. However, it is suboptimal to downsample randomly. Instead, in our first building block we want to downsample the point cloud while preserving contour, color and other feature information. This makes computation in the following blocks substantially cheaper since their complexities are polynomial in the number of points. We plan on parallelizing the graph-based approach described here. 
After downsampling our point cloud, our next building block segments this point cloud. There are multiple approaches to do this, but in our first implementation, we want to use the mean-shift or the quick-shift algorithm to perform the segmentation. This approach is again parallelizable. 

### Challenges
 In our CUDA approach, we plan on partitioning our volume into 3d cubes/voxels such that each thread block corresponds to
 a voxel. Because of the nature of the computation, there needs to be sharing of information between neighboring voxels. Minimizing inefficiencies here is a major challenge in maximizing speed. 
 
  Workload Imbalance: There is workload imbalance because of the uneven density of points. In CUDA, this corresponds to block-level imbalance. In addition, in the segmentation block, there is also some thread-level imbalance because of conditional statements. Minimizing execution divergence because of this is also a challenge.
  In both the segmentation and sampling blocks, the computation done for each point is spatially local. We can take advantage of this by designing efficient data structures and memory access patterns.
  
  At the same time, while tradeing off speed, we still want to preserve accuracy in our segmentation.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/rohanvarma16/pcseg/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
