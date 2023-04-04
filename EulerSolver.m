function [u,Kbc] = EulerSolver(nodes,conn,h,BC,P)
n = size(nodes,1);
K = zeros(n*6);
for i = 1:size(conn,1)
    x1 = nodes(conn(i,1),1);
    y1 = nodes(conn(i,1),2);
    x2 = nodes(conn(i,2),1);
    y2 = nodes(conn(i,2),2);
    L = sqrt((x2-x1)^2+(y2-y1)^2);
    alpha = atan2(y2-y1,x2-x1);
    Ki = BeamF(L,h(i),alpha);
    ki1 = 6*(conn(i,1)-1)+1;
    ki2 = 6*(conn(i,2)-1)+1;
    K([ki1:ki1+5,ki2:ki2+5],[ki1:ki1+5,ki2:ki2+5]) = K([ki1:ki1+5,ki2:ki2+5],[ki1:ki1+5,ki2:ki2+5]) + Ki;
end

idx = 1:6*n;
idx(BC) = [];
P = P(idx);
Kbc = K(idx,idx);

u = zeros(6*n,1);
us = Kbc\P;
u(idx,1) = us;
u = transpose(reshape(u,[6,n]));
end