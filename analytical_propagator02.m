function [tspan, oe_vec, ss_vec] =  analytical_propagator02(oe0, tspan, mu)
% This function computes the analytical evolution of the orbital elements 
% under the J2 approximation, following the Gauss perturbations eqns.
% 
% INPUTS
%   - oe0    : Vector of initial orbital elements ordered as 
%       - oe0(1) : a     - semi-major axis              [m]
%       - oe0(2) : e     - orbit eccentricity           [-]
%       - oe0(3) : i     - inclination                  [rad]
%       - oe0(4) : omega - argument of perigee          [rad]
%       - oe0(5) : Omega - RAAN                         [rad]
%       - oe0(6) : M     - mean anomaly                 [rad]
%   - tspan  : Vector of time properties                [s]
%   - mu     : Central body gravitational parameter     [m^3/s^2]
% 
% OUTPUTS
%   - tspan  : Vector of time properties                 [s]
%   - ss_vec : Final cartesian coordinates vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - oe_vec : Final keplerian coordinates vector (1x6) ordered as 
%       - oe0(1) : a     - semi-major axis              [m]
%       - oe0(2) : e     - orbit eccentricity           [-]
%       - oe0(3) : i     - inclination                  [deg]
%       - oe0(4) : w     - argument of perigee          [deg]
%       - oe0(5) : W     - RAAN                         [deg]
%       - oe0(6) : theta - true anomaly                 [deg]
% 
% REFERENCES
% Orbital Mechanics - Conway, Bruce A. and Prussing, John E. - 1993
% Chapter 7 - Numerical and analytical propagation - Pr. G. Kerschen
% see http://www.s3l.be/en/education

% Constants of the problem
R = 6371900;           % Earth's average radius (UGGI)     [m]
J = 1.082629e-3;       % Adimensional J2 term              [-]

% Extraction of initial keplerian elements
a   = oe0(1);
e   = oe0(2); 
i   = oe0(3);
w   = oe0(4);
W   = oe0(5);
M   = oe0(6);

% Initialisation of the coordinates vector
ss_vec = zeros(length(tspan), 6);
oe_vec = zeros(length(tspan), 6);
oe_vec(1,:) = oe0;

theta = M + ...
    (2 * e - e^3/4) * sin(M)   ...
    + 5/4   * e^2 * sin(2 * M) ...
    + 13/12 * e^3 * sin(3 * M) ;

oe_vec(1,end) = theta;

% Time increment
dt = tspan(2) - tspan(1);

Mdot = sqrt(mu/a^3);

for ii = 2:length(tspan)
    
    n = sqrt(mu/a^3);
    
    % Mean anomaly
    % Eq sl. 32 of Chapter 7 - Numerical and analytical propagation
    M = wrapTo2Pi(M + dt * Mdot);
    

    % Eccentric anomaly (use of Newton-Raphson)
    E = M; 
    diff = 1;
    while diff > 1e-5
        newE = wrapTo2Pi( ...
            E - (E - e * sin(E) - M) / ( 1 - e * cos(E) ) ...
                        );
        diff = abs(E - newE);
        E = newE;
    end
    
    
    % True anomaly
    theta = 2 * atan2(               ...
        sqrt(1 + e) * sin(E/2)   ,   ...
        sqrt(1 - e) * cos(E/2)   );

    % True anomaly trigonometry functions
    ct = cos(theta);
    st = sin(theta);
    
    % Central body distance
    r = a * ( 1 - e * cos(E) );
    
    % Disturbing acceleration
    
    F(1) = ( 1 - 3  * sin(i)^2 * st^2 ) / 2 ;
    F(2) = sin(i)^2 * st       * ct         ; 
    F(3) = sin(i)   * cos(i)   * st         ;
    
    F = - 3 * mu * J * R^2 / r^4 * F;
    
    
    % Extensively used term
    term = sqrt( a * ( 1 - e^2) / mu );
    
    % Eq 9.18
    adot = 2 * sqrt( a^3 / mu / (1 - e^2) ) * ( ...
        F(1) * e * st + F(2) * ( 1 + e * ct )          );
    
    % Eq. 9.22
    edot = term * ( F(1) * st + F(2) * ( ct + cos(E) ) );
    
    % Eq 9.28
    idot = term * F(3) * cos( theta ) / ( 1 + e * ct ) ;
        
    % Eq 9.36
    Wdot = term * F(3) * sin( theta ) / sin(i) / ( 1 + e * ct ) ;
    
    % Semilatus rectum
    p = a * (1 - e^2);

    Mdot = n * (1 + 3 / 4 * J * R^2 / p^2 * sqrt(1-e^2) * ( 2 - 3 * sin(i)^2));
    
    wdot = 3/4*J*R^2/p^2*(4 - 5 * sin(i)^2) * Mdot;
    
    a = a + dt * adot;
    e = e + dt * edot;
    i = i + dt * idot;
    w = w + dt * wdot;
    W = W + dt * Wdot;
    
    oe_vec (ii,1) = a;
    oe_vec (ii,2) = e;
    oe_vec (ii,3) = i;
    oe_vec (ii,4) = w;
    oe_vec (ii,5) = W;
    oe_vec (ii,6) = theta;
    
    % Transformation into cartesian coordinates
    ss_vec(ii,:) = kepl2cart_DF(oe_vec(ii,:), mu);
end
    
    % Transformation of angles units into degrees
    oe_vec(:,3:6) = rad2deg(oe_vec(:,3:6));

    
end