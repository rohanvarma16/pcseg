function G = graph_calcimp(G_inp, ptCloud)
    G = G_inp;
    for i=1:ptCloud.Count
       disp(i);
       xi_vec = [G.Nodes.X(i),G.Nodes.Y(i),G.Nodes.Z(i),...
           G.Nodes.R(i),G.Nodes.G(i),G.Nodes.B(i)]';
       neighbors_i = neighbors(G,int2str(i));
       sum_xj = zeros(length(xi_vec),1);
       for j=1:length(neighbors_i)
           j_ind = str2num(neighbors_i{j});
           xj_vec = [G.Nodes.X(j_ind),G.Nodes.Y(j_ind),G.Nodes.Z(j_ind),...
           G.Nodes.R(j_ind),G.Nodes.G(j_ind),G.Nodes.B(j_ind)]';
           weight = G.Edges.Weight(findedge(G,int2str(i),neighbors_i{j}));
           sum_xj = sum_xj + (weight .* xj_vec);
       end
    imp_i = norm(xi_vec - sum_xj)^2;
    G.Nodes.Imp(i) = imp_i;
    end
    
    %normalize importance:
    G.Nodes.Imp = G.Nodes.Imp/sum(G.Nodes.Imp);
end