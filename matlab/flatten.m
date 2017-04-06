map =parent_x;
for i=1:length(map)
   if(map(i)==0)
       map(i) = i;
   end
end
while 1
  map_ = map(map) ;
  if isequal(map_,map) ; break ; end
  map = map_ ;
end

[drop,drop,C] = unique(map)  ;

counts = zeros(max(C),1);
for i = 1:max(C)
counts(i) = sum(C==i);
end


[~,cluster_sort_idx] = sort(counts,'descend');

xyz_seg=[];
color_seg = uint8([]);
for i = 1:length(counts)
    %disp(i);
    cluster_ind = cluster_sort_idx(i);
    cluster_idx = find(C == cluster_ind);
    %new point cloud to show segment:
    xyz_seg_i =  [G_rs.Nodes.X(cluster_idx),G_rs.Nodes.Y(cluster_idx),G_rs.Nodes.Z(cluster_idx)];
    xyz_seg = [xyz_seg; xyz_seg_i];
    
    color_seg_i = 255* [G_rs.Nodes.R(cluster_idx),G_rs.Nodes.G(cluster_idx),G_rs.Nodes.B(cluster_idx)];
    color_avg_i = sum(color_seg_i)/length(cluster_idx);
    color_seg = uint8([color_seg;repmat(color_avg_i,length(cluster_idx),1)]);
    ptCloud_seg = pointCloud(xyz_seg, 'Color',color_seg);
    pcshow(ptCloud_seg);
    display_line = ['Cluster index: ',num2str(i),'Cluster size: ', num2str(length(cluster_idx))];
    disp(display_line);
     pause;
end


