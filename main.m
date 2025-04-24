%% Unknownmap Localcom 3D V1
% 

%%
clc;
clear;

%% 
%for TIME=1:50
%name = ['./case study robustness/map2 v1.2 robust/map2_', num2str(TIME), '.mat']; 
%load(name);
close all;
%% Network parameter

% Monitor area
%Obstacle_Area = genarea();
load("Obstacle_Area.mat");
%%
Covered_Area = zeros(size(Obstacle_Area,1),size(Obstacle_Area,2),size(Obstacle_Area,3));
[obs_x, obs_y, obs_z] = ind2sub(size(Obstacle_Area),find(Obstacle_Area==1));

% nodes info
MaxIt = 200;              % Maximum Number of Iterations
a = 1;                    % Acceleration Coefficient Upper Bound
N = 60;
rc = 30;
rs = 15;
sink=[90 90 90];
trap_thresh = 10;         % trap node condition
float_thresh= 100;        % float node condition

%moving parameter
v=5;                      % max velocity of node

%% Init first pop
figure;
initpop=unifrnd(90,100,[N 3]);
initpop(:,1:2)=unifrnd(70,100,[N 2]);
initpop(1,:)=sink;
pop=initpop;

% Node role Counter
no_move_counts = zeros(N,1);
no_profit_move_counts = zeros(N,1);
trap_matrix = zeros(N,1);

% Array to Hold Best Cost Values
BestCostIt = zeros(MaxIt, 1);
popIt=zeros(MaxIt,3*N);

%%     ABC Main Loop
for it = 1:MaxIt
    G=Graph(pop,rc);
    %% Trap matrix update
    al_trap_matrix=zeros(N,1);
    for l=2:N
        K=neighbors(G,l);
        al_trap_matrix(l) = no_move_counts(l) + mean(trap_matrix(K));
    end
    trap_matrix=al_trap_matrix;
    clear al_trap_matrix l K;

    %% Decision order
    % order of nodes decide which one will take move first

    distance = distances(G, 1,'Method','unweighted');
    N_layers = max(distance); % max layers= max distance to start node
    
    nodesInLayers = cell(N_layers, 2);  % batch layer matrix
    orderInLayers = cell(N_layers, 1);  % order layer matrix

    % separate to batchs
    for d = 1 : N_layers
        nodesInLayer = find(distance == d);
        sub_G=subgraph(G,nodesInLayer);
        sub_batchs= conncomp(sub_G);
        for i = 1: max(sub_batchs)
            k = sub_batchs==i;
            nodesInLayers{d,i} = nodesInLayer(k);  % Lưu các node vào cell tương ứng với từng tầng
        end
    end
    % Reorder
    for i = 1 : size(nodesInLayers,1)
        empty_order=[];
        maxLength = max(cellfun(@numel, nodesInLayers(i,:)));
        for k = 1 : maxLength
            for j = 1 : size(nodesInLayers,2)
                if numel(nodesInLayers{i,j}) >= k
                    empty_order=[empty_order nodesInLayers{i,j}(k)];
                end
            end
        end
        orderInLayers{i} = empty_order;
    end
    clear distance N_layers d sub_batchs sub_G i j k nodesInLayer maxLength empty_order nodesInLayers;
    
    %% Decision making
    for Layers = 1 : size(orderInLayers,1)
    for decision = 1:numel(orderInLayers{Layers})

        % node i makes move decision this turn
        i = orderInLayers{Layers}(decision);           
        al_pop=pop(i,:);             % alternative array of pop to load new positions
        G=Graph(pop,rc);

        % neighbor sensor group of node i
        K = neighbors(G,i);
        neighbor_pop=[];
        for n=1:numel(K)
            neighbor_pop=[neighbor_pop; pop(K(n),:)];
        end

        %% New Positions are created depends on types of node
        while (1)
            phi = a*unifrnd(-1, +1, [1 3])*(1-no_profit_move_counts(i)/MaxIt)^2;
            % --------------------------This node is trap node-----------------
            if no_move_counts(i) > trap_thresh
                [new_neigh_Cov, ~]=Cov_Func([neighbor_pop ]        ,rs,Obstacle_Area,Covered_Area);
                [old_neigh_Cov, ~]=Cov_Func([neighbor_pop;pop(i,:)],rs,Obstacle_Area,Covered_Area);
                % trap node can explore more
                if old_neigh_Cov-new_neigh_Cov>0.05*old_neigh_Cov
                    fitness_ratio=1;
                    %k = K(randi([1 numel(K)]));
                    al_pop = pop(i,:) + phi.*([v v v]);
                    no_profit_move_counts(i)=float_thresh/5;
                % trap node get stuck in local optimum
                else
                    fitness_ratio=0;
                    al_pop = pop(i,:) + phi.*([v v v]);
                    no_profit_move_counts(i)=no_profit_move_counts(i)+float_thresh/5;
                    no_move_counts(i)=0;
                end
            % --------------------------This node is float node------------------
            elseif no_profit_move_counts(i) > float_thresh
                if rand() >= 0.2
                    [~, n] = max(trap_matrix(K));
                    k=K(n);                             % k is the node that have the most value of no move counts
                    diff=abs(phi).*( pop(k,:) - pop(i,:) );
                    diff=max(min(diff,v),-v);
                    al_pop = pop(i,:) + diff;
                    fitness_ratio=0.993;
                else
                    al_pop = pop(i,:) + phi.*[v v v];
                    fitness_ratio=1;
                end
            % --------------------------This node is default node---------------
            else
                fitness_ratio=1;
                k = K(randi([1 numel(K)]));
                diff=(phi).*( pop(k,:) - pop(i,:) );
                diff=max(min(diff,v),-v);
                al_pop = pop(i,:) + diff;
            end
        %% boundary check
            al_pop(1) = min(max(al_pop(1,1), min(obs_x)+1),size(Obstacle_Area,1));
            al_pop(2) = min(max(al_pop(1,2), min(obs_y)+1),size(Obstacle_Area,2));
            al_pop(3) = min(max(al_pop(1,3), min(obs_z)+1),size(Obstacle_Area,3));
            obs_check1=[obs_x obs_y obs_z]-round(al_pop);
            obs_check2=abs(obs_check1(:,1))+abs(obs_check1(:,2))+abs(obs_check1(:,3));
            if ~any(obs_check2==0)
                break;
            end

        end
        %% Comparision of cost function
        neighbor_G=Graph([neighbor_pop; al_pop],rc);
        if Connectivity_graph(neighbor_G,[])==1
            [new_neigh_Cov, ~]=Cov_Func([neighbor_pop; al_pop  ],rs,Obstacle_Area,Covered_Area);
            [old_neigh_Cov, ~]=Cov_Func([neighbor_pop; pop(i,:)],rs,Obstacle_Area,Covered_Area);
            if (new_neigh_Cov) > (old_neigh_Cov)
                no_move_counts(i)=0;
                no_profit_move_counts(i) = max(no_profit_move_counts(i)-1,0);
                pop(i,:) = al_pop;
            elseif (new_neigh_Cov) > (old_neigh_Cov*fitness_ratio)
                no_move_counts(i)=0;
                no_profit_move_counts(i)= no_profit_move_counts(i)+1;
                pop(i,:) = al_pop;
            else
                no_move_counts(i) = no_move_counts(i)+1;
            end
        else
            no_move_counts(i) = no_move_counts(i)+1;
        end

    end
    end
    clear i j k K n al_pop neighbor_pop neighbor_G phi decisions_order decision new_neigh_Cov old_neigh_Cov fitness_ratio obs_check2 obs_check1;
    clear orderInLayers Layers diff;
    
    % Store Best Cost in that iteration
    [BestCostIt(it), ~] = Cov_Func(pop,rs,Obstacle_Area,Covered_Area);
    popIt(it,:)=reshape(pop,[1 N*3]);
    %disp([num2str(BestCostIt(it)) '  at iteration:  '  num2str(it)]);

    %% plot
    clf();
    hold on;
    plot3(pop(:,2),pop(:,1),pop(:,3),'ro','MarkerSize', 3,'Color','red')
    
    [x1,y1,z1] = sphere;
    for i=1:size(pop,1)
        x=x1*rs;
        y=y1*rs;
        z=z1*rs;
        surf(y+pop(i,2),x+pop(i,1),z+pop(i,3),'LineStyle',':','EdgeColor','cyan','EdgeAlpha',0.6,'FaceColor','none');
    end
    
    axis([0 100 0 100 0 100]); % Set the limits for X, Y, and Z axes
    
    [obs_x, obs_y, obs_z] = ind2sub(size(Obstacle_Area),find(Obstacle_Area==1));
    plot3(obs_y, obs_x, obs_z,'.', 'MarkerSize', 2, 'Color', 'blue');
    isosurface(0:100, 0:100, 0:100, Obstacle_Area, 0.5); % Correct dimension matching
    axis equal;
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(['3D Terrain Iteration: ' num2str(it) ', Coverage: ' num2str(BestCostIt(it))]);
    view(3);
    grid on;
    % Add lighting for better visualization
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting gouraud;
    drawnow;
    clear x y z x1 y1 z1 i ;
end
clear obs_x obs_y obs_z;
%%
%save(name)
%end
clear a ans Covered_Area float_thresh G initpop MaxIt no_move_counts no_profit_move_counts popIt2 trap_matrix trap_thresh;