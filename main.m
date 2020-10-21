%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AERO0024 - ASTRODYNAMICS %%%
%%%%%% ULiège - 2020-2021 %%%%%%

%%% Authors
% Axel DECHAMPS - S164160
% Casimir  FAYT - S196244


%%%%%% Orbital propagator %%%%%%
% This project consists in the development of an orbital propagator 
% of increasing complexity. The central body is the Earth.

close all; clear all; clc; format long;

%% Constants
mu = 398600.4418e9;       % Earth gravitational parameter [m^3/s^2]     

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

%% Initial ISS orbital parameters
e_ISS = .001379;                    % Eccentricity      [-]
a_ISS = 6794.57e3;                  % Semi-major axis   [m]
i_ISSd = 51.6445;                   % Inclination       [deg]
i_ISSr = deg2rad(i_ISSd);           % Inclination       [rad]
RAAN_ISSd = 128.8777;               % RAAN              [deg]
RAAN_ISSr =  deg2rad(RAAN_ISSd);    % RAAN              [rad]
omega_ISSd = 25.5173;               % Perigee argument  [deg]
omega_ISSr = deg2rad(omega_ISSd);   % Perigee argument  [rad]
M_ISSd = 146.2321;                  % Mean anomaly      [deg]
M_ISSr = deg2rad(M_ISSd);           % Mean anomaly      [rad]


% True anomaly from mean anomaly
% Fourier expansion from wikipedia
% See https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
theta_ISSr = M_ISSr + ...
    (2*e_ISS - e_ISS^3/4) * sin(M_ISSr) ...
    + 5/4   * e_ISS^2 * sin(2 * M_ISSr)   ...
    + 13/12 * e_ISS^3 * sin(3 * M_ISSr);

theta_ISSd = rad2deg(theta_ISSr);

oe_ISSd = [a_ISS, e_ISS, i_ISSd, omega_ISSd, RAAN_ISSd, theta_ISSd];
oe_ISSr = [a_ISS, e_ISS, i_ISSr, omega_ISSr, RAAN_ISSr, theta_ISSr];

% SL3 orbital propagator
[~, oe_SL3, ~, ce_SL3] = orbprop(oe_ISSd, 'time', tmax, 'dt', dt, 'fmodel', [0 0 0 0 0]);

%% Two-body propagator %%

% Iterations following Kepler equation
[~, oe_KEPL, ce_KEPL] = propagator01_KEPL_DECHAMPS_FAYT(oe_ISSr, tspan, mu);

% Numerical integration of Kepler relative motion
[~, oe_ODE, ce_ODE]  =  propagator01_ODE_DECHAMPS_FAYT(oe_ISSr, tspan, mu);


cartesian_comparison(ce_KEPL, ce_ODE, ce_SL3, tspan, MATLABc);

keplerian_comparison(oe_KEPL, oe_ODE, oe_SL3, tspan, MATLABc);


%% J2-term (Earth's oblateness) %%


%% Earth's atmosphere %%


%% Comparison with actual satellite data %%








%% Other functions
function cartesian_comparison(vec_KEPL, vec_ODE, vec_SL3, tspan, MATLABc)

f = figure;
f.WindowState = 'maximized';

subplot(3,2,1)
plot( tspan/3600 , vec_KEPL(:,1)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,1)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,1)/1000 , 'Color' , MATLABc{3}); 
title('x'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,3)
plot( tspan/3600 , vec_KEPL(:,2)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,2)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,2)/1000 , 'Color' , MATLABc{3}); 
title('y'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,5)
plot( tspan/3600 , vec_KEPL(:,3)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,3)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,3)/1000 , 'Color' , MATLABc{3}); 
title('z'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,2);
plot( tspan/3600 , vec_KEPL(:,4)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,4)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,4)/1000 , 'Color' , MATLABc{3}); 
title('xdot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');

subplot(3,2,4);
plot( tspan/3600 , vec_KEPL(:,5)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,5)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,5)/1000 , 'Color' , MATLABc{3}); 
title('ydot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');

subplot(3,2,6);
plot( tspan/3600 , vec_KEPL(:,6)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,6)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,6)/1000 , 'Color' , MATLABc{3});
title('zdot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');


leg = legend('Kepler propagator', 'ode45 integrator', 'SL3 propagator');
set(leg,'Position', [.455 .01 .125 .075],'Units', 'normalized');

end


function keplerian_comparison(vec_KEPL, vec_ODE, vec_SL3, tspan, MATLABc)

f = figure;
f.WindowState = 'maximized';

subplot(3,2,1);
plot( tspan/3600 , vec_KEPL(:,1)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,1)/1000 , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,1)/1000 , 'Color' , MATLABc{3}); 
title('Semi-major axis'); ylabel('a [km]'); xlabel('Time [hours]');

subplot(3,2,3);
plot( tspan/3600 , vec_KEPL(:,2) , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_ODE(:,2) , 'Color' , MATLABc{2});
plot( tspan/3600 ,  vec_SL3(:,2) , 'Color' , MATLABc{3}); 
title('Eccentricity'); ylabel('e [-]'); xlabel('Time [hours]');

subplot(3,2,5);
plot( tspan/3600 , rad2deg(vec_KEPL(:,3)) , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  rad2deg(vec_ODE(:,3)) , 'Color' , MATLABc{2});
plot( tspan/3600 ,           vec_SL3(:,3) , 'Color' , MATLABc{3}); 
title('Inclination'); ylabel('i [deg]'); xlabel('Time [hours]');

subplot(3,2,2);
plot( tspan/3600 , rad2deg(vec_KEPL(:,4)) , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  rad2deg(vec_ODE(:,4)) , 'Color' , MATLABc{2});
plot( tspan/3600 ,           vec_SL3(:,4) , 'Color' , MATLABc{3}); 
title('Argument of perigee'); 
ylabel('\omega [deg]'); xlabel('Time [hours]');

subplot(3,2,4);
plot( tspan/3600 , rad2deg(vec_KEPL(:,5)) , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  rad2deg(vec_ODE(:,5)) , 'Color' , MATLABc{2});
plot( tspan/3600 ,           vec_SL3(:,5) , 'Color' , MATLABc{3}); 
title('RAAN'); ylabel('\Omega [deg]'); xlabel('Time [hours]');

subplot(3,2,6);
plot( tspan/3600 , rad2deg(vec_KEPL(:,6)) , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  rad2deg(vec_ODE(:,6)) , 'Color' , MATLABc{2}); 
plot( tspan/3600 ,           vec_SL3(:,6) , 'Color' , MATLABc{3});
title('True anomaly'); ylabel('\theta [deg]'); xlabel('Time [hours]');


leg = legend('Kepler propagator', 'ode45 integrator', 'SL3 propagator');
set(leg,'Position', [.455 .01 .125 .075],'Units', 'normalized');

end


