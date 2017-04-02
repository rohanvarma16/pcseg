function G = graph_addnodes(ptCloud)

nrows = size(ptCloud.Location,1);
T = table();
G = graph();
parfor i=1:nrows
    %disp(i);
    nodeprops = {{int2str(i)},ptCloud.Location(i,1),...
    ptCloud.Location(i,2),ptCloud.Location(i,3),...
    ptCloud.Color(i,1),ptCloud.Color(i,2), ...
    ptCloud.Color(i,3),'VariableNames', {'NodeID',...
    'x','y','z','r','g','b'}};
    T = [T ; nodeprops];
end
    G = addnode(G, T);         
end


       
