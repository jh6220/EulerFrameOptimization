function [u,ubt,loadb,Pe,Me,Se] = EulerBucklingSolver(nodes,conn,h,BC,P)
[u,K] = EulerSolver(nodes,conn,h,BC,P);
[Pe,Me,Se] = PostProcessingEuler3D(nodes,conn,u,h);
n = size(nodes,1);
Kgeo = zeros(n*6);
for i = 1:size(conn,1)
    x1 = nodes(conn(i,1),1);
    y1 = nodes(conn(i,1),2);
    x2 = nodes(conn(i,2),1);
    y2 = nodes(conn(i,2),2);
    L = sqrt((x2-x1)^2+(y2-y1)^2);
    alpha = atan2(y2-y1,x2-x1);
    Kgeoi = EulerBeamGeoK(L,alpha,Pe(i));
    ki1 = 6*(conn(i,1)-1)+1;
    ki2 = 6*(conn(i,2)-1)+1;
    Kgeo([ki1:ki1+5,ki2:ki2+5],[ki1:ki1+5,ki2:ki2+5]) = Kgeo([ki1:ki1+5,ki2:ki2+5],[ki1:ki1+5,ki2:ki2+5]) + Kgeoi;
end

idx = 1:6*n;
idx(BC) = [];
Kgeo = Kgeo(idx,idx);
[ubs,loadb] = eig(K,-Kgeo);
ub = zeros(6*n,size(ubs,2));
ub(idx,:) = ubs;
loadb = diag(loadb);
ubt = zeros(n,6,size(ub,2));
for i = 1:size(ub,2)
ubt(:,:,i) = transpose(reshape(ub(:,i),[6,n]));
end
[~,idx] = sort(abs(loadb));
loadb = loadb(idx);
ubt = ubt(:,:,idx);

% for i = BC
%     Kgeo(i,:)=0;Kgeo(:,i)=0;
%     Kgeo(i,i)=1;
%     K(i,:)=0;K(:,i)=0;
%     K(i,i)=1;
% end
% [ub,loadb] = eig(K,Kgeo);
% loadb = diag(loadb);
% ubt = zeros(n,6,size(ub,2));
% for i = 1:size(ub,2)
% ubt(:,:,i) = transpose(reshape(ub(:,i),[6,n]));
% end
end