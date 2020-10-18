function [oe_vec] = cart2kepl_KZ(ss_vec, mu)
% This function transforms state-space vector into oe_vecerian elements

% INPUTS
%   - ss_vec : State-space vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - mu     : Central body gravitational parameter     [km^3/s^2]

% OUTPUTS
%   - oe_vec : Vector of Keplerian elements (1x6) ordered as follows :
%       - oe_vec(1) : a     - semi-major axis       [m]
%       - oe_vec(2) : e     - orbit eccentricity    [-]
%       - oe_vec(3) : i     - inclination           [rad]
%       - oe_vec(4) : omega - argument of perigee   [rad]
%       - oe_vec(5) : Omega - RAAN                  [rad]
%       - oe_vec(6) : theta - true anomaly          [rad]

% REFERENCES
% Chapter 3A - The orbit in space of Pr. G. Kerschen
% see http://www.s3l.be/en/education

I = [1 0 0];    % Towards vernal equinoxe
K = [0 0 1];    % Normal to equatorial plane
J = [0 1 0];    % Perpendicular to I and K

N = length(ss_vec(1,:));

oe_vec = zeros(6,N);

for ii = 1:N

    %% Position and velocity vectors
    r = ss_vec(1:3,ii);    rnorm = norm(r);
    v = ss_vec(4:end,ii);  vnorm = norm(v);

    %% Specific angular momentum
    h = cross(r,v);     hnorm = norm(h);

    %% Eccentricity
    e = norm(cross(v,h) - mu * r / rnorm) / mu;

    %% Semi-major axis
    a = rnorm / ( 2 - rnorm * vnorm^2 / mu);

    %% Inclination
    i = acos( h(3) / hnorm );

    %% Longitude
    n = cross( K, h / hnorm );  % Nodal vector
    nnorm = norm(n);

    % Vector with direction corresponding to apse line
    apse = cross(v,h) / mu - r / rnorm;
    apsenorm = norm(apse);

    RAAN = acos ( dot( n, I ) / nnorm );
    if dot(n,J) < 0
        RAAN = 2 * pi - RAAN;
    end

    %% Argument of perigee
    omega = acos ( dot(n, apse)  / nnorm / apsenorm );
    if dot(apse,K) < 0
        omega = 2 * pi - omega;
    end

    %% True anomaly
    theta = acos ( dot(r,apse) / rnorm / apsenorm );
    if dot(r,v) < 0
        theta = 2 * pi - theta;
    end

    oe_vec(:,ii) = [
        a;
        e;
        i;
        omega;
        RAAN;
        theta];
    
end
end