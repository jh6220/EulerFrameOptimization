function [] = PlotStructure(nodes,conn)
figure(1)
scatter(nodes(:,1),nodes(:,2))
for i = 1:size(conn,1)
    hold on
    plot([nodes(conn(i,1),1),nodes(conn(i,2),1)],[nodes(conn(i,1),2),nodes(conn(i,2),2)],'k')
end
axis equal
end