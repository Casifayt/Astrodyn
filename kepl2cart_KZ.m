function [ss_vec] = kepl2cart_KZ(oe, mu)
% This function transforms keplerian elements into cartesian elements
%
% INPUTS
%   - oe_vec : Vector of Keplerian elements (1x6) ordered as follows :
%       - oe_vec(1) : a     - semi-major axis           [m]
%       - oe_vec(2) : e     - orbit eccentricity        [-]
%       - oe_vec(3) : i     - inclination               [rad]
%       - oe_vec(4) : omega - argument of perigee       [rad]
%       - oe_vec(5) : Omega - RAAN                      [rad]
%       - oe_vec(6) : theta - true anomaly              [rad]
%   - mu     : Earth gravitational parameter            [m^3/s^2]
%
% OUTPUTS
%   - ss_vec : State-space vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
% 
% REFERENCES
% From Rene Schwarz materials
% See https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf


% Extraction of initial keplerian elements
a     = oe(1);
e     = oe(2);
i     = oe(3);
omega = oe(4);
RAAN  = oe(5);
theta = oe(6);

% Eccentric anomaly
E = 2 * atan2 ( sqrt(1-e) * sin(theta/2), sqrt(1+e) * cos(theta/2) );

% Mean anomaly
M = E - e * sin(E);

% Central body distance
r = a * ( 1 - e * cos(E) );

% Position and velocity vector in orbital frame
pos_orb = r * [
    cos(theta);
    sin(theta);
    0];

vel_orb = sqrt(mu * a) / r * [
    -sin(E);
    sqrt(1-e^2) * cos(E);
    0];

% Transformation into inertial frame
% Rotation matrices
Rz_RAAN = [
       cos(-RAAN)  -sin(-RAAN)    0     ;
       sin(-RAAN)   cos(-RAAN)    0     ;
            0           0         1     ];

Rx_i = [
            1           0         0     ;
            0       cos(-i)    -sin(-i) ;     
            0       sin(-i)     cos(-i) ];

Rz_per = [
       cos(-omega)  -sin(-omega)    0     ;
       sin(-omega)   cos(-omega)    0     ;
            0           0           1     ];

% Cartesian position in inertial frame
r_cart    = Rz_RAAN * Rx_i * Rz_per * pos_orb;

% Cartesian velocity in inertial frame
rdot_cart = Rz_RAAN * Rx_i * Rz_per * vel_orb;

% State-space vector
ss_vec = [r_cart; rdot_cart];

end







