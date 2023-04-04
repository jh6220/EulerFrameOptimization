function [P_d,M_b,S] = PostProcessingEuler3D(nodes,conn,u,h)
E = 70000;
A = 6.35.*h;
Iz = 6.35*h.^3/12;

n_e = size(conn,1);
P_d = zeros(n_e,1);
S = zeros(n_e,1);
M_b = zeros(n_e,2);

for i = 1:n_e
    Ui = transpose(cat(2,u(conn(i,1),:),u(conn(i,2),:)));
    x1 = nodes(conn(i,1),1);
    y1 = nodes(conn(i,1),2);
    x2 = nodes(conn(i,2),1);
    y2 = nodes(conn(i,2),2);
    L = sqrt((x2-x1)^2+(y2-y1)^2);
    alpha = atan2(y2-y1,x2-x1);
    
    c = cos(alpha);
    s = sin(alpha);

    T = [c, s,  0,  0,  0,  0;
        -s, c,  0,  0,  0,  0;
        0,  0,  1,  0,  0,  0;
        0,  0,  0,  1,  0,  0;
        0,  0,  0,  0,  1,  0;
        0,  0,  0,  0,  0,  1];
    
    T = [T,zeros(6,6);zeros(6,6),T];
    
    ui = T*Ui;
    P_d(i) = E*A(i)*(ui(7)-ui(1))/L;
    M_b(i,1) = E*Iz(i)*(6*(ui(8)-ui(2))-L*(4*ui(6)+2*ui(12)))/L^2;
    M_b(i,2) = E*Iz(i)*(6*(-ui(8)+ui(2))+L*(2*ui(6)+4*ui(12)))/L^2;
    S(i) = E*Iz(i)*(12*(ui(2)-ui(8))+6*L*(ui(6)+ui(12)))/L^3;
end
end