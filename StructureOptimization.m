% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         30,-100];
% 
% conn = [1,  2;
%         1,  3;
%         1,  4;
%         2,  4;
%         3,  4];

% nodes = [0         0
%  -140.0000   30.0000
%   320.0000   30.0000
%    85.9942  111.5616];

% h = [20,20,20,20,20];

% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         0,-100;
%         120,-100];
% 
% conn = [1,  2;
%         1,  3;
%         1,  4;
%         1,  5;
%         2,  4;
%         4,  5;
%         5,  3];
% 
% h = 20*ones(1,size(conn,1));
% 
% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         95,   31;
%         0,  101;
%         104, 93];

% 
% nodes = [0, 0;
%         -140, 30;
%         320,  30;
%         120,   20;
%         0,  101;
%         10, 100];

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
% 
nodes = [0, 0;
        -140, 30;
        320,  30;
        0,    -60;
        0,  100];

conn = [2,  4;
        2,  5;
        4,  1;
        5,  1;
        3,  4;
        3,  5];

h = 20*ones(1,size(conn,1));

PlotStructure(nodes,conn)
%%
nodes_best = nodes;
h_best = h;

BC = [8,9,10,11,13,14,15,16,17];
P = zeros(size(nodes,1),6);
P(1,2) = -8000;
P = reshape(transpose(P),[6*size(nodes,1),1]);
dL = 50;
w = [10^3,10^0];
ds = -10^-5;
d_step = [10^-4.5,10^-4];
n_steps = 5000;
Cost_history = zeros(n_steps,1);
log_arr = zeros(n_steps,3);

figure(1)
anim = animatedline;
[Cost,~] = StructureCost(nodes,conn,h,BC,P,dL,w);
Cost_min = Cost;
addpoints(anim,1,Cost);
for i_step = 1:n_steps
    [Cost,log] = StructureCost(nodes,conn,h,BC,P,dL,w);
    Cost_history(i_step) = Cost;
    log_arr(i_step,:) = log;
    if (mod(i_step,1) == 0) && i_step>1
        disp(log)
%         disp(nodes(4:end,:))
        addpoints(anim,i_step,Cost_history(i_step));
        axis([0 i_step min(Cost_history(1:i_step)) max(Cost_history(1:i_step))])
        drawnow
    end

    if Cost<Cost_min
        nodes_best = nodes;
        h_best = h;
        Cost_min = Cost;
    elseif (Cost>1.1*Cost_min)
        disp("break")
        break
    end
    
    Cost_nodes = zeros(size(nodes,1)-3,2);
    for i = 1:(size(nodes,1)-3)
        for j = 1:2
            d_nodes = zeros(size(nodes));
            d_nodes(i+3,j) = ds;
            [Cost_nodes(i,j),~] = StructureCost(nodes+d_nodes,conn,h,BC,P,dL,w);
        end
    end

    Nabla_nodes = (Cost_nodes-Cost)/ds*d_step(2);
    Nabla_nodes(Nabla_nodes>1) = 1;
    Nabla_nodes(Nabla_nodes<-1) = -1;
    Nabla_nodes = [zeros(3,2);Nabla_nodes];

    nodes = nodes - Nabla_nodes;

    for i = 1:5
        Nabla_h = Nabla_hF(nodes,conn,h,BC,P,dL,w,ds)*d_step(1);
        Nabla_h(Nabla_h>0.1) = 0.1;
        Nabla_h(Nabla_h<-0.1) = -0.1;
        h = h - Nabla_h;
        h(h<3.5) = 3.5;
    end
end

%%
PlotStructure(nodes,conn)
title("Surface Area: "+num2str(log(1)*1000)+"mm^2")

%%
i = '6';
writematrix(h,['Structures/h_',i,'.csv'])
writematrix(nodes,['Structures/nodes_',i,'.csv'])
writematrix(conn,['Structures/conn_',i,'.csv'])

%%
i = '6';
h = readmatrix(['Structures/h_',i,'.csv']);
nodes = readmatrix(['Structures/nodes_',i,'.csv']);
conn = readmatrix(['Structures/conn_',i,'.csv']);