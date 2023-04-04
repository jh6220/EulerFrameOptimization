function [Cost,log] = StructureCost(nodes,conn,h,BC,P,dL,w)

stress_y = 260-2;

[u,~,loadb,~,~,Pe,Me,~,conn_full_idx,L] = EulerSolverMeshing(nodes,conn,h,BC,P,dL);
Me_max = zeros(size(conn,1),1);
Pe_max = zeros(size(conn,1),1);
for i = 1:size(conn,1)
    Me_max(i) =  max(abs(Me(conn_full_idx(i,1):conn_full_idx(i,2),:)),[],'all');
    Pe_max(i) =  max(abs(Pe(conn_full_idx(i,1):conn_full_idx(i,2),:)),[],'all');
end
stress_max = max(Pe_max./(6.35.*h') + 6*Me_max./(h'.^2.*6.35));
% h_min = transpose((Pe_max+sqrt(Pe_max.^2+24*Me_max*6.35*stress_y))/(2*6.35*stress_y));
% h = max([h_min;h;ones(size(h))*3]);
loadb_max = loadb(find(loadb>0,1));
% uy00 = u(1,2);
% Cost = h*L + 1000*w(1)^(1.5-loadb_max) + w(2)^(stress_max-stress_y);
Cost = h*L + (w(1)*max(1.005-loadb_max,0))^2 + (w(2)*max(stress_max-stress_y,0))^2;
% Cost = h*L + (w(1)*max(1.1-loadb_max,0)) + (w(2)*max(stress_max-stress_y,0));
log = [h*L/1000,loadb_max,stress_max];
% log = [h*L,(w(1)*max(1.1-loadb_max,0))^2 ,(w(2)*max(stress_max-stress_y,0))^2];
end