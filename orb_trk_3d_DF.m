function orb_trk_3d_DF (cart_vec, tspan, rot)
% This function plots the orbit around the Earth in 3D
%
% INPUTS
%   - cart_vec  : Array of cartesian coordinates        [m] 
%   - tspan     : Vector of time                        [s]
%   - rot       : Parameter to animate the 3D Earth     [-]
%
% Copyrights for the Earth 3D plot to Ryan Gray, 2013
% Taken from MathWorks at 
% https://nl.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example

alpha = .8;         % Transparency of the globe

% Plot the 3D Earth
earth3d(alpha);

% Plot the 3D orbit
x_vec = cart_vec(:,1);
y_vec = cart_vec(:,2);
z_vec = cart_vec(:,3);

plot3(x_vec, y_vec, z_vec, 'LineWidth', 2.5);

plot3(x_vec(1), y_vec(1), z_vec(1), 's', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(x_vec(1) + 1e5, y_vec(1) + 1e5, z_vec(1) + 1e4, ...
    'Start', 'Color', 'r');


plot3(x_vec(end), y_vec(end), z_vec(end), 'o', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(x_vec(end) + 1e5, y_vec(end) + 1e5, z_vec(end) + 1e3, ...
    'End', 'Color', 'r', 'HorizontalAlignment', 'right');

end

function earth3d(alpha)
%% Textured 3D Earth example
%
% Ryan Gray
% 8 Sep 2004
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013

%% Options

space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.

image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

% Mean spherical earth

erad    = 6371008.7714; % equatorial radius (meters)
prad    = 6371008.7714; % polar radius (meters)

%% Create figure

figure('Color', space_color, 'Name', '3D Orbit', 'WindowState', 'maximized');

hold on;

% Turn off the normal axes

set(gca,...
    'NextPlot','add',       ...
    'Visible','off',        ...
    'Interactions',         ...
    [rulerPanInteraction    ...
    zoomInteraction         ...
    rotateInteraction       ...
    dataTipInteraction      ...
    ]);

axis equal; axis auto;

% Set initial view
view(0,30);

axis vis3d;

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);


%% Texturemap the globe

% Load Earth image for texture map

cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');


end
