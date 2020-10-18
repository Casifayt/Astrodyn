function [ss_vec] = propagator01_KEPL_DECHAMPS_FAYT(oe0, tspan, mu)
% This function computes the evolution of the orbital elements by solving the
% Kepler equation, using the Newton-Raphson method

% Kepler equation
%       M(t) = E(t) - e * sin(E);

% Newton-Raphson method
%       E(j+1) = E(j) - f( E(j) ) / f'( E(j) )

% INPUTS
%   - oe0   : Vector of initial orbital elements ordered as follows
%       - oe0(1) : a     - semi-major axis          [m]
%       - oe0(2) : e     - orbit eccentricity       [-]
%       - oe0(3) : i     - inclination              [rad]
%       - oe0(4) : omega - argument of perigee      [rad]
%       - oe0(5) : Omega - RAAN                     [rad]
%       - oe0(6) : theta - true anomaly             [rad]
%   - tspan : Time vector                           [s]
%   - mu    : Central body gravitational parameter  [km^3/s^2]

% OUTPUTS
%   - ss_vec : State-space vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]

% REFERENCES
% From Rene Schwarz materials
% See https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

% Extraction of initial keplerian elements
a       = oe0(1);
e       = oe0(2);
i       = oe0(3);
omega   = oe0(4);
RAAN    = oe0(5);
theta0  = oe0(6);


% Initialisation of the orbital parameters vectors
M       = zeros(1,length(tspan));   % Mean anomaly              [rad]
theta   = zeros(1,length(tspan));   % True anomaly              [rad]
E       = zeros(1,length(tspan));   % Eccentricity anomaly      [rad]
r       = zeros(1,length(tspan));   % Distance                  [m]
pos_orb = zeros(3,length(tspan));   % Position in orbital frame [m]
vel_orb = zeros(3,length(tspan));   % Velocity in orbital frame [m/s]

% Initial conditions
M(1)    = theta0;
theta(1)= theta0;

% Time increment
dt = tspan(2) - tspan(1);

for ii = 2:length(tspan)
    % Mean anomaly
    M(ii) = wrapTo2Pi(M(ii-1) + dt * sqrt(mu/a^3));

    % Eccentric anomaly (use of Newton-Raphson)
    E(ii) = M(ii); 
    diff = 1;
    while diff > 1e-5
        newE = E(ii) - (E(ii) - sin(E(ii)) - M(ii)) / ( 1 - cos(E(ii)) );
        diff = abs(E(ii) - newE);
        E(ii) = newE;
    end

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
    
    
    
    
    
    
    
    