function [tspan, oe_vec, ss_vec] = propagator01_KEPL_DECHAMPS_FAYT(oe0, tspan, mu)
% This function computes the evolution of the orbital elements by solving the
% Kepler equation, using the Newton-Raphson method
% 
% C�cile equation
%       M(t) = E(t) - e * sin(E);
% 
% Newton-C�cile method
%       E(j+1) = E(j) - f( E(j) ) / f'( E(j) )
% 
% INPUTS
%   - oe0    : Vector of initial orbital elements ordered as follows
%       - oe0(1) : a     - semi-major axis              [m]
%       - oe0(2) : e     - orbit eccentricity           [-]
%       - oe0(3) : i     - inclination                  [rad]
%       - oe0(4) : omega - argument of perigee          [rad]
%       - oe0(5) : Omega - RAAN                         [rad]
%       - oe0(6) : theta - true anomaly                 [rad]
%   - tspan  : Vector of time properties                 [s]
%   - mu     : Central body gravitational parameter      [m^3/s^2]
% 
% OUTPUTS
%   - tspan  : Vector of time properties                 [s]
%   - ss_vec : Final cartesian coordinates vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - oe_vec : Final keplerien coordinates vector (1x6) ordered as the
%   input vector
% 
% 
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
% Eccentric anomaly
E(1) = 2 * atan2 ( sqrt(1-e) * sin(theta0/2), sqrt(1+e) * cos(theta0/2) );

% Mean anomaly
M(1) = E(1) - e * sin(E(1));

theta(1)= theta0;

% Time increment
dt = tspan(2) - tspan(1);

for ii = 1:length(tspan)
    if ii == 1
        oe_vec = oe0;
        ss_vec = kepl2cart_KZ(oe0, mu);

        oe_vec = oe_vec';
        ss_vec = ss_vec';
    else
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
        Rz_RAAN = [
               cos(RAAN)  -sin(RAAN)    0     ;
               sin(RAAN)   cos(RAAN)    0     ;
                    0           0         1   ];

        Rx_i = [
                    1           0         0    ;
                    0       cos(i)    -sin(i)  ;     
                    0       sin(i)     cos(i) ];

        Rz_per = [
               cos(omega)  -sin(omega)    0     ;
               sin(omega)   cos(omega)    0     ;
                    0           0           1  ];

        % Cartesian position in inertial frame
        r_cart    = Rz_RAAN * Rx_i * Rz_per * pos_orb;

        % Cartesian velocity in inertial frame
        rdot_cart = Rz_RAAN * Rx_i * Rz_per * vel_orb;

        % State-space vector
        ss_vec = [r_cart; rdot_cart];

        % Orbital elements vector
        oe_vec = cart2kepl_KZ(ss_vec, mu);
        ss_vec = ss_vec';
        oe_vec = oe_vec'; 
    end
end




end
    
    
    
    
    
    
    
    