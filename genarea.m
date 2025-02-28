function Obstacle_Area = genarea()
% Interest area value =0
% Obstacle area value =1
bound=100;
% Define grid size (X and Y coordinates)
x = 0:bound; % X-axis range
y = 0:bound; % Y-axis range
z = 0:bound; % Z-axis range (depth levels)
%[X, Y, Z] = meshgrid(x, y, z); % 3D grid for the entire volume

%% Generate random terrain (Z values) using Perlin-like noise
rng(1); % Seed for reproducibility
[X2D, ~] = meshgrid(x, y); % 2D grid for terrain
terrain = zeros(size(X2D)); % Initialize terrain

%% Create fractal noise for terrain
scales = [10, 20, 40]; % Noise scales for detail levels
amplitudes = [10, 5, 2]; % Heights for each noise layer
for i = 1:length(scales)
    % Generate random noise in the range [0, 1]
    noise = rand(ceil(size(X2D, 1)/scales(i)), ceil(size(X2D, 2)/scales(i)));
    
    % Interpolate this noise to fit the grid size of X2D
    noise = interp2(noise, linspace(1, size(noise, 1), size(X2D, 1)), ...
                    linspace(1, size(noise, 2), size(X2D, 2))');
    terrain = terrain + amplitudes(i) * noise; % Accumulate layers
end

%% Normalize terrain height to fit the range [0,50]
terrain = terrain - min(terrain(:)); % Set minimum height to 0
terrain = terrain / max(terrain(:)) * (bound/2); % Scale terrain height to [0.50]

%% Generate 3D bitmap
Obstacle_Area = zeros(length(y), length(x), length(z)); % Initialize 3D binary matrix
for k = 1:length(z)
    Obstacle_Area(:, :, k) = terrain >= z(k); % Fill voxels below the terrain
end
%{
% Visualize the 3D bitmap using isosurface
figure;
axis([0 100 0 100 0 100]); % Set the limits for X, Y, and Z axes
isosurface(X, Y, Z, Obstacle_Area, 0.5); % Correct dimension matching
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Terrain Spanning X-Y Plane, Depth [-100, -50]');
view(3);
grid on;

% Add lighting for better visualization
light('Position', [1 1 1], 'Style', 'infinite');
lighting gouraud;
%}


