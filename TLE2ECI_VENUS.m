function [oe_start,oe_final,ss_start,ss_final] = TLE2ECI_VENUS(mu)

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
%       - oe_vec(3) : i     - inclination           [rad]
%       - oe_vec(4) : omega - argument of perigee   [rad]
%       - oe_vec(5) : Omega - RAAN                  [rad]
%       - oe_vec(6) : theta - true anomaly          [rad]

%% Conversion of the TLEs (we only retain TLE n°1 and TLE n°4)

[r_vec,v_vec] = tle2eci('tle_Venus.txt'); 
ss_start = [r_vec(1,1) ; r_vec(1,2) ; r_vec(1,3) ; 
            v_vec(1,1) ; v_vec(1,2) ; v_vec(1,3)];
ss_final = [r_vec(4,1) ; r_vec(4,2) ; r_vec(4,3) ; 
            v_vec(4,1) ; v_vec(4,2) ; v_vec(4,3)];  

%% Conversion of the cartesian coordinates into the keplerian coordinates

oe_start = cart2kepl_KZ(ss_start,mu);
oe_start(3:6) = rad2deg(oe_start(3:6));
oe_final = cart2kepl_KZ(ss_final,mu);
oe_final(3:6) = rad2deg(oe_final(3:6));

end