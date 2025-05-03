function [coverage, Covered_Area] = Cov_Func_v2(pop,rs,Obstacle_Area,Covered_Area)
%This is fitness function to cal random area coverage ratio

%pop is a Nx3dim matrix holding position of nodes
%rs is the sensing rad of nodes
%Obstacle_Area is a 3dim matrix   
%Covered_Area =0 : uncovered area of interest
%Covered_Area =1 : covered area of interest
%Covered_Area =-2: covered area of obstacle
%Covered_Area =-1: uncovered area of obstacle


%% recover sensor uncovered area

%find all covered point and recover them to uncover status
[obs_x, obs_y, obs_z] = ind2sub(size(Covered_Area),find(Covered_Area==1));
for i = 1:numel(obs_x)
    Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = 0;
end

%% check sensor covered area
inside_sector = false(size(Covered_Area,1),size(Covered_Area,2),size(Covered_Area,3));
%%
for j=1:(size(pop,1))
    %%
    % Node position j-th
    x0 = pop(j,1);
    y0 = pop(j,2);
    z0 = pop(j,3);
    rsJ=rs(j);

    % Boundary constraint
    x_ub=min(max(ceil(x0+rsJ),1),size(Covered_Area,1));
    x_lb=min(max(ceil(x0-rsJ),1),size(Covered_Area,1));
    y_ub=min(max(ceil(y0+rsJ),1),size(Covered_Area,2));
    y_lb=min(max(ceil(y0-rsJ),1),size(Covered_Area,2));
    z_ub=min(max(ceil(z0+rsJ),1),size(Covered_Area,3));
    z_lb=min(max(ceil(z0-rsJ),1),size(Covered_Area,3));

    % Local Grid
    [X, Y, Z] = meshgrid(linspace(x_lb, x_ub, x_ub-x_lb+1), linspace(y_lb, y_ub, y_ub-y_lb+1), linspace(z_lb, z_ub, z_ub-z_lb+1));
    
    % Distance matrix
    D = sqrt((X - x0).^2 + (Y - y0).^2 + (Z - z0).^2);
    
    % In rs condition
    in_circle = D <= rsJ;
    
    % distance conditions
    inside_sector(y_lb:y_ub,x_lb:x_ub,z_lb:z_ub) = inside_sector(y_lb:y_ub,x_lb:x_ub,z_lb:z_ub) | (in_circle); 
    
end       
%%
Covered_Area = inside_sector.* (1-Obstacle_Area);
%clear D Theta in_circle in_angle inside_sector;


%% add obstacle to covered area
%% check obstacle in all covered area
[obs_x, obs_y, obs_z] = ind2sub(size(Obstacle_Area),find(Obstacle_Area==1));
for i = 1:numel(obs_x)
    if Covered_Area (obs_x(i), obs_y(i), obs_z(i)) == 1   
        Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = -2;
    end
end

count1=numel(ind2sub(size(Covered_Area),find(Covered_Area==1)));		                          % count covered points on wanted location  (wanted)
%count2=numel(ind2sub(size(Covered_Area),find(Covered_Area==-2)));		                          % count covered points on unwanted location (obstacles)
count3=numel(Obstacle_Area)-numel(ind2sub(size(Obstacle_Area),find(Obstacle_Area==1)));           % count total points on wanted location

%coverage=((count1-count2)/count3);	    % function to avoid obstacles
coverage=(count1/count3);		        % function to aim on wanted area

%% recover obs covered area
% 
% [obs_x, obs_y, obs_z] = ind2sub(size(Covered_Area),find(Covered_Area==-2));
% for i = 1:numel(obs_x)
%     Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = -1;
% end
