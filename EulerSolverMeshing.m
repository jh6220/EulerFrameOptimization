function [u,ub,loadb,nodes_full,conn_full,Pe,Me,Se,conn_full_idx,L,h_full] = EulerSolverMeshing(nodes,conn,h,BC,P,dL)
n = size(nodes,1);
nodes_full = nodes;
conn_full = [];
h_full = [];
conn_full_idx = zeros(size(conn));
L = zeros(size(conn,1),1);
for i = 1:size(conn,1)
    n1 = nodes(conn(i,1),:);
    n2 = nodes(conn(i,2),:);
    L(i) = sqrt((n1(1)-n2(1))^2+(n1(2)-n2(2))^2);
%     if h(i)>10
%         nn=6;
%     else
%         nn=2;
%     end
%     nn = max(ceil(L(i)/dL),2);
    nn = 4;
    if nn>2
        scaler = reshape(linspace(0,1,nn),nn,1);
        nodes_add = n1+(n2-n1).*scaler(2:end-1);
        nn_count = size(nodes_full,1);
        conn_add = zeros(nn-1,2);
        conn_add(1,:) = [conn(i,1),nn_count+1];
        for j = 2:nn-1
            conn_add(j,:) = [nn_count+j-1,nn_count+j];
        end
        conn_add(end,:) = [nn_count+nn-2,conn(i,2)];
    else
        conn_add = conn(i,:);
        nodes_add = [];
    end
    h_add = h(i)*ones(nn-1,1);
    nodes_full = cat(1,nodes_full,nodes_add);
    conn_full_idx(i,:) = [size(conn_full,1)+1,size(conn_full,1)+size(conn_add,1)];
    conn_full = cat(1,conn_full,conn_add);
    h_full = cat(1,h_full,h_add);
end
P_full = cat(1,P,zeros((size(nodes_full,1)-n)*6,1));
% if plotIf
%     PlotStructure(nodes_full,conn_full);
%     hold on
% end
% [u_full] = TimoshenkoSolver(nodes_full,conn_full,h_full,BC,P_full);
% 
% [s_d_full,s_b_full,s_xz_full] = PostProcessing(nodes_full,conn_full,u_full,h_full);

% if plotIf
%     PlotDeformed(nodes_full,conn_full,u_full,s_d_full)
% end

[u,ub,loadb,Pe,Me,Se] = EulerBucklingSolver(nodes_full,conn_full,h_full,BC,P_full);

end