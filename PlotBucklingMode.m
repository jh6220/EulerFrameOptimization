function [] = PlotBucklingMode(nodes,conn,ub,iMode)
figure(2)
nodes = [nodes,zeros(size(nodes,1),1)];
scatter3(nodes(:,3),nodes(:,1),nodes(:,2))
for i = 1:size(conn,1)
    hold on
    plot3([nodes(conn(i,1),3),nodes(conn(i,2),3)],[nodes(conn(i,1),1),nodes(conn(i,2),1)],[nodes(conn(i,1),2),nodes(conn(i,2),2)],'k')
end

nodes = nodes+ub(:,1:3,iMode)*50;
% nodes = nodes+u(:,1:3)/max(abs(u(:,1:3)))*10;
hold on
scatter3(nodes(:,3),nodes(:,1),nodes(:,2))
for i = 1:size(conn,1)
    hold on
    plot3([nodes(conn(i,1),3),nodes(conn(i,2),3)],[nodes(conn(i,1),1),nodes(conn(i,2),1)],[nodes(conn(i,1),2),nodes(conn(i,2),2)],'r')
end

axis equal
end