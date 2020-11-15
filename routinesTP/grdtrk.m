function grdtrk(rECEF, str, rot, dt)
% This function displays the ground track of the spacecraft.
% INPUTS
%   - rECEF : Array of Cartesian coordinates [Nx3]  [m]
%   - str   : Title of the plot                     [-]
%   - rot   : Rotation of the Earth (0 or 1)        [-]
%   - dt    : Time step of the simulation           [s] 
%
% Copyrights to the SL3 propagator for most of the script

earth_ang_vel = 360 / 86400;    % Earth angular vel [deg/s]
n = size(rECEF, 1);

% Initialisation of geodetic coordinates arrays
LON = zeros(n, 1);
LAT = zeros(n, 1);
H = zeros(n, 1);

% Use of given function to compute geodetic coordinates
for j = 1 : n
    [H(j), LON(j), LAT(j)] = ecef2geodetic(rECEF(j, :)');
end

% Conversion to degrees and meters
LAT = rad2deg(LAT);
LON = rad2deg(LON);
Hgdtrk = H * 1e-3;

% Computation of the Earth's rotation (if asked)
if rot == 1
    for j = 1:length(LON)
        if j > 0
            LON(j) = LON(j) + j * dt * earth_ang_vel;
        end
        if LON(j) > 180
            LON(j) = LON(j) - 360;
        end
    end
end

% Find discontinuites
N = length(LAT);
for j = 1 : (N - 1)
    if abs(LON(j + 1) - LON(j)) > 300
        LON = [LON(1 : j); NaN; LON((j + 1) : N)];
        LAT = [LAT(1 : j); NaN; LAT((j + 1) : N)];
        Hgdtrk = [Hgdtrk(1 : j); NaN; Hgdtrk((j + 1) : N)];
        N = N + 1;
    end
end

%% Plot
% Font size
set(0, 'defaultaxesfontsize', 16); set(0, 'defaulttextfontsize', 16);

% Display of plot
box on; axis on;

% Earth's surface
image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
data = imread(image_file);

image('CData', data, 'XData', [-180 180], 'YData', [90 -90]);
hold on; grid on;

% Groundtrack plot
plot(LON, LAT, 'r', 'linewidth', 0.2);

% Initial and final markers
in = plot(LON(1), LAT(1), 's', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
fin = plot(LON(end), LAT(end), 'o', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
legend([in, fin], 'Start', 'End');

% Plot title
title(['Ground track - ' str]);
% X axis
xlabel('Longitude [deg]'); xlim([-180 180]);
% Y axis
ylabel('Latitude [deg]'); ylim([-90 90]);
% Axis ticks
set(gca, 'xtick', -180 : 60 : 180, 'ytick', -90 : 30 : 90);

hold off;

return

function [h, lambda, phi] = ecef2geodetic(rECEF)

% 
% [h, lambda, phi] = ecef2geodetic(rECEF)
% 
% Provide the geodetic coordinates of a position vector in the ECEF frame.
% 
% Inputs:
%   - rECEF = 3-element column vector defining a position in the ECEF frame
%            [m]
% 
% Outputs:
%   - h = altitude from the reference ellipsoid [m]
%   - lambda = geodetic longitude [rad]
%   - phi = geodetic latitude [rad]
% 
% Ref: MATLAB function ecef2geodetic.
% 
% Lamberto Dell'Elce
% 

persistent a e2 ep2 f b

if isempty(a)
    % Ellipsoid constants
    a  = 6378.137e3; % Semimajor axis
    e2 = 0.081819190842^2; % Square of first eccentricity
    ep2 = e2 / (1 - e2); % Square of second eccentricity
    f = 1 - sqrt(1 - e2); % Flattening
    b = a * (1 - f); % Semiminor axis
end

x = rECEF(1);
y = rECEF(2);
z = rECEF(3);

%% Longitude
lambda = atan2(y,x); % [rad]

%% Latitude
% Distance from Z-axis
RHO = hypot(x,y);

% Bowring's formula for initial parametric (beta) and geodetic (phi)
% latitudes
beta = atan2(z, (1 - f) * RHO);
phi = atan2(z   + b * ep2 * sin(beta)^3,...
            RHO - a * e2  * cos(beta)^3);

% Fixed-point iteration with Bowring's formula
% (typically converges within two or three iterations)
betaNew = atan2((1 - f)*sin(phi), cos(phi));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    phi = atan2(z + b * ep2 * sin(beta)^3, ...
        RHO - a * e2  * cos(beta).^3); % [rad]
    betaNew = atan2((1 - f)*sin(phi), cos(phi));
    count = count + 1;
end

%% Altitude
% Calculate ellipsoidal height from the final value for latitude
sin_phi = sin(phi);
N = a / sqrt(1 - e2 * sin_phi^2);
h = RHO * cos(phi) + (z + e2 * N * sin_phi) * sin_phi - N; % [m]

return