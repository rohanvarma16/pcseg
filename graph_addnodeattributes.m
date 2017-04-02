function G_new = graph_addnodeattributes(G,ptCloud)
    G_new = G;
    G_new.Nodes.X = ptCloud.Location(:,1);
    G_new.Nodes.Y = ptCloud.Location(:,2);
    G_new.Nodes.Z = ptCloud.Location(:,3);
    G_new.Nodes.R = ptCloud.Color(:,1);
    G_new.Nodes.G = ptCloud.Color(:,2);
    G_new.Nodes.B = ptCloud.Color(:,3);
    G_new.Nodes.Imp = zeros(ptCloud.Count,1);
end