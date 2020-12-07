function [oe_SGP4,ss_SGP4] = SGP4_VENUS(mu,tspan)

% This function computes osculating position and velocity from a file of 
% TLEs. It then return the cartesian elements (ss) as well as the keplerian
% elements (oe) of the satellite.
% 
% INPUTS
%   - mu     : Central body gravitational parameter     [m^3/s^2]
%
% OUTPUTS
%   - ss_vec : State-space vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - oe_vec : Vector of Keplerian elements (1x6) ordered as follows :
%       - oe_vec(1) : a     - semi-major axis       [m]
%       - oe_vec(2) : e     - orbit eccentricity    [-]
%       - oe_vec(3) : i     - inclination           [deg]
%       - oe_vec(4) : omega - argument of perigee   [deg]
%       - oe_vec(5) : Omega - RAAN                  [deg]
%       - oe_vec(6) : theta - true anomaly          [deg]

%% Simulation of the TLEs with the SGP4 Propagator (we only retain TLE 1)

[r_vec,v_vec] = sgp4(tspan','tle_Venus.txt',1);
r_vec = r_vec';
v_vec = v_vec';
ss_SGP4 = [r_vec ; v_vec]';

%% Conversion of the cartesian coordinates into the keplerian coordinates

oe_SGP4 = cart2kepl_DF(ss_SGP4',mu)';
oe_SGP4(:,3:6) = oe_SGP4(:,3:6);

end