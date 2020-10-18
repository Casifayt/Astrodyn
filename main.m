%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AERO0024 - ASTRODYNAMICS %%%
%%%%%% ULiège - 2020-2021 %%%%%%

%%% Authors
% Axel DECHAMPS - S164160
% Casimir  FAYT - S196244


%%%%%% Orbital propagator %%%%%%
% This project consists in the development of an orbital propagator 
% of increasing complexity. The gravitational body is the Earth.

close all; clear all; clc;

%% Constants
mu = 398600.4418e9;       % Earth gravitational parameter [m^3/s^2]     

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
[~, oe_SL3, ~, ce_SL3] = orbprop(oe_ISSd, 'time', tmax, 'dt', dt);

%% Two-body propagator %%

[~, oe_KEPL, ce_KEPL] = propagator01_KEPL_DECHAMPS_FAYT(oe_ISSr, tspan, mu);

[~, oe_ODE, ce_ODE]  =  propagator01_ODE_DECHAMPS_FAYT(oe_ISSr, tspan, mu);

figure;
subplot(3,2,1)
plot(tspan,ce_KEPL(:,1),'r'); hold on;
plot(tspan,ce_ODE(:,1),'b');
plot(tspan,ce_SL3(:,1),'k'); title('x');

subplot(3,2,3)
plot(tspan,ce_KEPL(:,2),'r'); hold on;
plot(tspan,ce_ODE(:,2),'b');
plot(tspan,ce_SL3(:,2),'k'); title('y');

subplot(3,2,5)
plot(tspan,ce_KEPL(:,3),'r'); hold on;
plot(tspan,ce_ODE(:,3),'b');
plot(tspan,ce_SL3(:,3),'k'); title('z');

subplot(3,2,2);
plot(tspan,ce_KEPL(:,4),'r'); hold on;
plot(tspan,ce_ODE(:,4),'b');
plot(tspan,ce_SL3(:,4),'k'); title('xdot');

subplot(3,2,4);
plot(tspan,ce_KEPL(:,5),'r'); hold on;
plot(tspan,ce_ODE(:,5),'b');
plot(tspan,ce_SL3(:,5),'k'); title('ydot');

subplot(3,2,6);
plot(tspan,ce_KEPL(:,6),'r'); hold on;
plot(tspan,ce_ODE(:,6),'b');
plot(tspan,ce_SL3(:,6),'k'); title('zdot');




%% J2-term(Earth's oblateness) %%


%% Earth's atmosphere %%


%% Comparison with actual satellite data %%





