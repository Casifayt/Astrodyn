function [tspan, oe_vec, ss_vec] =  analytical_propagator01(oe0, tspan, mu)
% This function computes the analytical evolution of the orbital elements 
% by solving the Kepler equation, using the Newton-Raphson method
% 
% Kepler equation
%       M(t) = E(t) - e * sin(E);
% 
% Newton-Raphson method
%       E(j+1) = E(j) - f( E(j) ) / f'( E(j) )
% 
% INPUTS
%   - oe0    : Vector of initial orbital elements ordered as follows
%       - oe0(1) : a     - semi-major axis              [m]
%       - oe0(2) : e     - orbit eccentricity           [-]
%       - oe0(3) : i     - inclination                  [rad]
%       - oe0(4) : omega - argument of perigee          [rad]
%       - oe0(5) : Omega - RAAN                         [rad]
%       - oe0(6) : M     - mean anomaly                 [rad]
%   - tspan  : Vector of time properties                 [s]
%   - mu     : Central body gravitational parameter      [m^3/s^2]
% 
% OUTPUTS
%   - tspan  : Vector of time properties                 [s]
%   - ss_vec : Final cartesian coordinates vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - oe_vec : Final keplerian coordinates vector (1x6) ordered as the
%   input vector
% 
% 
% REFERENCES
% From Rene Schwarz materials
% See https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

% Extraction of initial keplerian elements
a   = oe0(1);
e   = oe0(2);
i   = oe0(3);
w   = oe0(4);
W   = oe0(5);
M0  = oe0(6);


% Initialisation of the orbital parameters vectors
M       = zeros(1,length(tspan));   % Mean anomaly              [rad]
theta   = zeros(1,length(tspan));   % True anomaly              [rad]
E       = zeros(1,length(tspan));   % Eccentricity anomaly      [rad]
r       = zeros(1,length(tspan));   % Distance                  [m]
pos_orb = zeros(3,length(tspan));   % Position in orbital frame [m]
vel_orb = zeros(3,length(tspan));   % Velocity in orbital frame [m/s]

M(1) = M0;

% Time increment
dt = tspan(2) - tspan(1);

for ii = 1:length(tspan)
    if ii ~= 1
        % Mean anomaly
        M(ii) = wrapTo2Pi(M(ii-1) + dt * sqrt(mu/a^3));
    end

    % Eccentric anomaly (use of Newton-Raphson)
    E(ii) = M(ii); 
    diff = 1;
    while diff > 1e-5
        newE = E(ii) - (E(ii) - e * sin(E(ii)) - M(ii)) / ( 1 - e * cos(E(ii)) );
        diff = abs(E(ii) - newE);
        E(ii) = newE;
    end
    
    E(ii) = wrapTo2Pi(E(ii));

    % True anomaly
    theta(ii) = 2 * atan2(               ...
        sqrt(1+e) * sin(E(ii)/2)     ,   ...
        sqrt(1-e) * cos(E(ii)/2)     );

    % Central body distance
    r(ii) = a * ( 1 - e * cos(E(ii)) );

    % Position and velocity vectors in orbital frame
    pos_orb(:,ii) = r(ii) * [
        cos(theta(ii));
        sin(theta(ii));
        0];

    vel_orb(:,ii) = sqrt(mu * a) / r(ii) * [
        -sin(E(ii));
        sqrt(1-e^2) * cos(E(ii));
        0];

    % Transformation into inertial frame
    % Rotation matrices
    Rz_W = [
           cos(W)  -sin(W)    0     ;
           sin(W)   cos(W)    0     ;
                0         0     1     ];

    Rx_i = [
                1           0         0    ;
                0       cos(i)    -sin(i)  ;     
                0       sin(i)     cos(i) ];

    Rz_w = [
           cos(w)  -sin(w)    0 ;
           sin(w)   cos(w)    0 ;
                0          0    1 ];

    % Cartesian position in inertial frame
    r_cart    = Rz_W * Rx_i * Rz_w * pos_orb;

    % Cartesian velocity in inertial frame
    rdot_cart = Rz_W * Rx_i * Rz_w * vel_orb;

    % State-space vector
    ss_vec = [r_cart; rdot_cart];

    % Orbital elements vector
    oe_vec = cart2kepl_DF(ss_vec, mu)';
    ss_vec = ss_vec';
end
end
    