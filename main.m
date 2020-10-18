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
tspan = dt:dt:tmax;

%% Initial ISS orbital parameters
e_ISS = .001379;                % Eccentricity      [-]
a_ISS = 6794.57e3;              % Semi-major axis   [m]
i_ISS = deg2rad(51.6445);       % Inclination       [rad]
RAAN_ISS = deg2rad(128.8777);   % RAAN              [rad]
omega_ISS = deg2rad(25.5173);   % Perigee arg       [rad]
M_ISS = deg2rad(146.2321);      % Mean anomaly      [rad]

oe_ISS = [a_ISS, e_ISS, i_ISS, omega_ISS, RAAN_ISS, M_ISS];

% SL3 orbital propagator
[t, oe_SL3] = orbprop([a_ISS; e_ISS; i_ISS; omega_ISS; RAAN_ISS; M_ISS], ...
    'time', tmax, ...
    'dt', dt);

%% Two-body propagator %%

ss_vec = propagator01_KEPL_DECHAMPS_FAYT(oe_ISS, tspan, mu);
oe_vec = cart2kepl_KZ(ss_vec(:,end), mu);

% ss_vec = propagator01_ODE_DECHAMPS_FAYT(oe_ISS, tspan, mu);
% oe_vec = cart2kepl(ss_vec(:,end), mu);

% plot(tspan,ss_vec,'-'); legend;

%% J2-term(Earth's oblateness) %%


%% Earth's atmosphere %%


%% Comparison with actual satellite data %%