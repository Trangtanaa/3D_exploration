function [coverage , Covered_Area] = Cov_Func_v1(pop,rs,Obstacle_Area,Covered_Area)
%This is fitness function to cal random area coverage ratio

%pop is a Nx3dim matrix holding position of nodes
%rs is the sensing rad of nodes
%Obstacle_Area is a 3dim matrix   
%Covered_Area =0 : uncovered area of interest
%Covered_Area =1 : covered area of interest
%Covered_Area =-2: covered area of obstacle
%Covered_Area =-1: uncovered area of obstacle


%% recover sensor uncovered area

% find all covered point and recover them to uncover status
[obs_x, obs_y, obs_z] = ind2sub(size(Covered_Area),find(Covered_Area==1));
for i = 1:numel(obs_x)
    Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = 0;
end

%% check all covered area
for j=1:size(pop,1)
    % x value of sensor j is pop(j,1)
    % y value of sensor j is pop(j,2)
    % z value of sensor j is pop(j,3)
    start_point=[floor(pop(j,1)) floor(pop(j,2)) floor(pop(j,3))];
    for i = 0:(rs+1)
        for k = 0:(rs+1)
            for l = 0:(rs+1)
                map_x= start_point(1);
                map_y= start_point(2);
                map_z= start_point(3);
                % center point : [mapx mapy mapz]
                dist = sqrt((map_x+i-pop(j,1))^2+(map_y+k-pop(j,2))^2+(map_z+l-pop(j,3))^2);
                % case dist <= rs
                if (dist < rs || dist== rs)
                    if map_x+i <= size(Obstacle_Area,1) && map_y+k <= size(Obstacle_Area,2) && map_z+l <=size(Obstacle_Area,3)
                        Covered_Area(map_x+i,map_y+k,map_z+l) = 1;
                    end
                    if map_x+i <= size(Obstacle_Area,1) && map_y+k <= size(Obstacle_Area,2) && map_z-l >0
                        Covered_Area(map_x+i,map_y+k,map_z-l) = 1;
                    end
                    if map_y+k <= size(Obstacle_Area,1) && map_x-i >0 && map_z+l <=size(Obstacle_Area,3)  
                        Covered_Area(map_x-i,map_y+k,map_z+l) = 1;
                    end
                    if map_y+k <= size(Obstacle_Area,1) && map_x-i >0 && map_z-l>0 
                        Covered_Area(map_x-i,map_y+k,map_z-l) = 1;
                    end
                    if map_y-k > 0 && map_x+i <= size(Obstacle_Area,2) && map_z+l <=size(Obstacle_Area,3)  
                        Covered_Area(map_x+i,map_y-k,map_z+l) = 1;
                    end
                    if map_y-k > 0 && map_x+i <= size(Obstacle_Area,2) && map_z-l>0  
                        Covered_Area(map_x+i,map_y-k,map_z-l) = 1;
                    end
                    if map_y-k > 0 && map_x-i > 0 && map_z+l <=size(Obstacle_Area,3)      
                        Covered_Area(map_x-i,map_y-k,map_z+l) = 1;
                    end
                    if map_y-k > 0 && map_x-i > 0 && map_z-l>0    
                        Covered_Area(map_x-i,map_y-k,map_z-l) = 1;
                    end
                end
            end
        end
    end
end

%% check obstacle in all covered area
[obs_x, obs_y, obs_z] = ind2sub(size(Obstacle_Area),find(Obstacle_Area==1));
for i = 1:numel(obs_x)
    if Covered_Area (obs_x(i), obs_y(i), obs_z(i)) == 1   
        Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = -2;
    end
end

count1=numel(ind2sub(size(Covered_Area),find(Covered_Area==1)));		                          % count covered points on wanted location  (wanted)
count2=numel(ind2sub(size(Covered_Area),find(Covered_Area==-2)));		                          % count covered points on unwanted location (obstacles)
count3=numel(Obstacle_Area)-numel(ind2sub(size(Obstacle_Area),find(Obstacle_Area==1)));           % count total points on wanted location

%coverage=((count1-count2)/count3);	    % function to avoid obstacles
coverage=(count1/count3);		        % function to aim on wanted area

%% recover obs covered area
%{
% find all obs points that are covered and recover them to uncover status
[obs_x, obs_y, obs_z] = ind2sub(size(Covered_Area),find(Covered_Area==-2));
for i = 1:numel(obs_x)
    Covered_Area(obs_x(i), obs_y(i), obs_z(i)) = -1;
end
%}
