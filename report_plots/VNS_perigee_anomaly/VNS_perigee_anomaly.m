%% Perigee and true anomaly erratic behaviours explanation
% This script plots the perigee and true anomaly of Venµs according to
% different propagators. These behaviours are of particular interest and
% needs special plot to explain it, as the true anomaly sometimes decrease,
% which is counterintuitive as one could conclude that the satellite goes
% backwards on its trajectory. It is linked to the quasi-curclarity of the
% orbit, yielding discontinuous jumps in the perigee argument, and as the
% true anomaly is defined in function of it, it can lead to "impossible"
% behaviour

% Credits to Michelle Hirsch for the ds2nfu function
% https://nl.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion

clear; close all; clc;

MATLABc = {
    [0, 0.4470, 0.7410];
    [0.8500, 0.3250, 0.0980];
    [0.9290, 0.6940, 0.1250];
    [0.4940, 0.1840, 0.5560];
    [0.4660, 0.6740, 0.1880];
    [0.3010, 0.7450, 0.9330];
    [0.6350, 0.0780, 0.1840];
    };

% Time properties
tmax = 86400;
N = 1e3;
dt = tmax/N;
tspan = 0:dt:tmax;

load DF3_1d_1e5.mat;
load S3L_1d_1e5.mat;
load SGP4_1d_1e5.mat;

% Derivative treshold
der_tres = 20;

% Indices for values to be plotted (focus on a few periods, as the graph is
% not clear when plotting full 24h propagation
start = 1;      % Initial index
fin = 75;       % Final index (approximately two hours of propagation
d = 1;          % Increment

% Reduction of time vector
tspan = tspan(start:d:fin)/3600;

%% From S3L propagator
% True anomaly (2pi phase shifts unwrapped)
f_S3Ld =  rad2deg(unwrap(deg2rad(oe_S3L (start:d:fin,6))));
% Numerical derivative (to identify fast varying zones)
f_S3L_deriv = gradient(f_S3Ld);

% Argument of perigee (2pi phase shifts unwrapped)
w_S3Ld =  rad2deg(unwrap(deg2rad(oe_S3L (start:d:fin,4))));
% Numerical derivatives
w_S3L_deriv = gradient(w_S3Ld);

%% From SGP4 propagator
% True anomaly (2pi phase shifts unwrapped)
f_SGP4d = rad2deg(unwrap(deg2rad(oe_SGP4(start:d:fin,6))));

% Numerical derivatives
f_SGP4_deriv = gradient(f_SGP4d);

% Argument of perigee (2pi phase shifts unwrapped)
w_SGP4d = rad2deg(unwrap(deg2rad(oe_SGP4(start:d:fin,4))));

% Numerical derivatives
w_SGP4_deriv = gradient(w_SGP4d);

%% From DF3 propagator
% True anomaly (2pi phase shifts unwrapped)
f_DF3d =  rad2deg(unwrap(deg2rad(oe_DF3 (start:d:fin,6))));

% Numerical derivatives
f_DF3_deriv = gradient(f_DF3d);

% Argument of perigee (2pi phase shifts unwrapped)
w_DF3d =  rad2deg(unwrap(deg2rad(oe_DF3 (start:d:fin,4))));

% Numerical derivatives
w_DF3_deriv = gradient(w_DF3d);


%% Plot
f = figure;

% Plot of the true anomaly evolution
subplot(2,1,1);
plot(tspan, f_S3Ld, 'LineWidth', 2); hold on;
plot(tspan, f_DF3d, 'LineWidth', 2);
plot(tspan, f_SGP4d, 'LineWidth', 2); box on; grid on;
ylabel('\theta [deg]');

% Identification of the inversion of sign of the derivative of the true
% anomaly evolution. The corresponding points coordinates in subplot axis
% will be converted to figure coordinates in order to plot annotating line

DF3_in = [];    % Initialisation of array of initial coordinates (DF3)
idx_DF3 = [];   % Array containing indices of values mathcing criteria (DF3)

SGP4_in = [];    % Initialisation of array of initial coordinates (SGP4)
idx_SGP4 = [];   % Array containing indices of values mathcing criteria (SGP4)

% If the sign of the previous derivative was positive and the sign of
% the current derivative is negative, we have peculiar point
for i = 1:length(tspan)
   if i ~= 1
        if sign(f_DF3_deriv(i-1)) == 1 && sign(f_DF3_deriv(i)) == -1
          % Checking for derivative values matching criteria for DF3 vector

          % Transformation from axis to figure coordinates
          [x_fin, y_fin] = ds2nfu(tspan(i-2), f_DF3d(i-2));

          % Storing corresponding index
          idx_DF3 = [idx_DF3 i-2];

          % Storing initial point of annotation in figure coordinates
          DF3_in = [DF3_in; x_fin, y_fin];
        end
   end
   
   if abs(f_SGP4_deriv(i)) > der_tres
      % Checking for derivative values matching criteria for SGP4 vector 
       
      % Transformation from axis to figure coordinates
      [x_fin, y_fin] = ds2nfu(tspan(i), f_SGP4d(i));
      
      % Storing corresponding index
      idx_SGP4 = [idx_SGP4 i];
      
      % Storing initial point of annotation in figure coordinates
      SGP4_in = [SGP4_in; x_fin, y_fin];
   end
end

% Plot of the argument of perigee evolution
subplot(2,1,2);

plot(tspan, w_S3Ld, 'LineWidth', 2); hold on;
plot(tspan, w_DF3d, 'LineWidth', 2);
plot(tspan, w_SGP4d, 'LineWidth', 2); box on; grid on;
ylabel('\omega [deg]');

% Identification of the points of interest in argument of perigee 
% evolution according to the index of the values showing derivative 
% inversion in the true anomaly plot

DF3_fin = [];  % Initialisation of array of final coordinates (DF3)

for i = 1:length(idx_DF3)
    % Transformation of axis into figure coordinates
    [x_in, y_in] = ds2nfu(tspan(idx_DF3(i)), w_DF3d(idx_DF3(i)));      
    
    % Storing final point of annotation in figure coordinates
    DF3_fin = [DF3_fin; x_in, y_in];
end

SGP4_fin = []; % Initialisation of array of final coordinates (SGP4)

for i = 1:length(idx_SGP4)
    % Transformation of axis into figure coordinates
    [x_in, y_in] = ds2nfu(tspan(idx_SGP4(i)), w_SGP4d(idx_SGP4(i)));      
    
    % Storing final point of annotation in figure coordinates
    SGP4_fin = [SGP4_fin; x_in, y_in];
end

% Plot of all annotations for DF3/S3L propagators
for i = 1:length(idx_DF3)
    annotation('line',          ...     % Line type
    [DF3_in(i,1) DF3_fin(i,1)], ...     % Initial and final x coordinates
    [DF3_in(i,2) DF3_fin(i,2)], ...     % Initial and final y coordinates
    'Color', 'r', 'LineStyle', '--'  ); % Red dashed line
end

% Plot of all annotations for SGP4 propagator
for i = 1:length(idx_SGP4)
    annotation('line',          ...     % Line type
    [SGP4_in(i,1) SGP4_fin(i,1)], ...     % Initial and final x coordinates
    [SGP4_in(i,2) SGP4_fin(i,2)], ...     % Initial and final y coordinates
    'Color', 'b', 'LineStyle', '--'  ); % Red dashed line
end

% Legend placement
leg = legend('S3L','DF3', 'SGP4');
set(leg,'Position', [.8 .48 .125 .075], 'Units', 'normalized');

    
