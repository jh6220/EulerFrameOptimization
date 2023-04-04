function [] = PlotDeformed(nodes,conn,u,c)
figure(1)
c = (c-min(c))/(max(c)-min(c));

max_x = max(nodes(:,1));
min_x = min(nodes(:,1));
max_y = max(nodes(:,2));
min_y = min(nodes(:,2));
L = max(max_x-min_x,max_y-min_y);
u_nodes = u(:,1:2)./max(abs(u(:,1:2)),[],'all')*0.1*L;
nodes = nodes+u_nodes;
axis equal
scatter(nodes(:,1),nodes(:,2))
for i = 1:size(conn,1)
    hold on
    plot([nodes(conn(i,1),1),nodes(conn(i,2),1)],[nodes(conn(i,1),2),nodes(conn(i,2),2)],'color',[c(i),0,1-c(i)])
end
end