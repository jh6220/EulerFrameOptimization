function [Nabla_h] = Nabla_hF(nodes,conn,h,BC,P,dL,w,ds)
Cost = StructureCost(nodes,conn,h,BC,P,dL,w);
Cost_h = zeros(size(h));
for i = 1:length(h)
    delta_h = zeros(size(h));
    delta_h(i) = ds;
    [Cost_h(i),~] = StructureCost(nodes,conn,h+delta_h,BC,P,dL,w);
end
Nabla_h = (Cost_h-Cost)/(ds);

end