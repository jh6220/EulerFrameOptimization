% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         56,135];
% 
% conn = [1,  2;
%         1,  3;
%         1,  4;
%         2,  4;
%         3,  4];
% 
% h = [20,20,20,20,20];

nodes = [0, 0;
        -140, 30;
        320,  30;
        9,98;
        136,99];

conn = [1,  2;
        1,  3;
        1,  4;
        1,  5;
        2,  4;
        4,  5;
        5,  3];

h = 2*[2.51,1.5,1.5,1.5,9.7,12.3,9.6];

% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         100,    0;
%         -30,  100;
%         120, 100];
% 
% conn = [1,  2;
%         1,  4;
%         4,  3;
%         2,  5;
%         5,  6;
%         6,  3;
%         1,  5;
%         4,  6];
% 
% h = 20*ones(1,size(conn,1));

%%

BC = [8,9,10,11,13,14,15,16,17];
P = zeros(size(nodes,1),6);
P(1,2) = -8000;
P = reshape(transpose(P),[6*size(nodes,1),1]);
dL = 20;
[u,ub,loadb,nodes_full,conn_full,Pe,Me,Se,conn_full_idx,L,h_full] = EulerSolverMeshing(nodes,conn,h,BC,P,dL);

Me_max = zeros(size(conn,1),1);
Pe_max = zeros(size(conn,1),1);
for i = 1:size(conn,1)
    Me_max(i) =  max(abs(Me(conn_full_idx(i,1):conn_full_idx(i,2),:)),[],'all');
    Pe_max(i) =  max(abs(Pe(conn_full_idx(i,1):conn_full_idx(i,2),:)),[],'all');
end
% stress_max = max(Pe_max./(6.35.*h') + 6*Me_max./(h'.^2.*6.35));
stress = (Pe_max./(6.35.*h') + 6*Me_max./(h'.^2.*6.35));
PlotStructure(nodes_full,conn_full)
hold on
PlotDeformed(nodes_full,conn_full,u,-Pe)
% PlotDeformed(nodes,conn,u(1:size(nodes),:),Pe_max./(h'*6.35))
i_b = find(loadb>0,1);
hold off
PlotBucklingMode(nodes_full,conn_full,ub,i_b)

disp(loadb(i_b))

%%
nodes = nodes_full;
conn = conn_full;
h = h_full';