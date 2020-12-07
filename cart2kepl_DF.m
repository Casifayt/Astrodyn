function [oe_vec] = cart2kepl_DF(ss_vec, mu)
% This function transforms state-space vector into oe_vecerian elements
%
% INPUTS
%   - ss_vec : State-space vector (1x6) ordered as follows :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]
%   - mu     : Central body gravitational parameter     [m^3/s^2]
%
% OUTPUTS
%   - oe_vec : Vector of Keplerian elements (1x6) ordered as follows :
%       - oe_vec(1) : a     - semi-major axis       [m]
%       - oe_vec(2) : e     - orbit eccentricity    [-]
%       - oe_vec(3) : i     - inclination           [deg]
%       - oe_vec(4) : omega - argument of perigee   [deg]
%       - oe_vec(5) : Omega - RAAN                  [deg]
%       - oe_vec(6) : theta - true anomaly          [deg]
%
% REFERENCES
% Chapter 3A - The orbit in space - Pr. G. Kerschen
% see http://www.s3l.be/en/education

I = [1 0 0];    % Towards vernal equinoxe
K = [0 0 1];    % Normal to equatorial plane
J = [0 1 0];    % Perpendicular to I and K

N = length(ss_vec(1,:));

oe_vec = zeros(6,N);

for ii = 1:N

    %% Position and velocity vectors
    r = ss_vec(1:3,ii);     rnorm = norm(r);
    v = ss_vec(4:end,ii);   vnorm = norm(v);

    %% Specific angular momentum
    h = cross(r,v);         hnorm = norm(h);

    %% Eccentricity
    e = norm(cross(v,h) - mu * r / rnorm) / mu;
    
    %% Semi-major axis
    a = rnorm / ( 2 - rnorm * vnorm^2 / mu);

    %% Inclination
    i = acos( h(3) / hnorm );

    %% Longitude
    n = cross( K, h / hnorm );  % Nodal vector
    nnorm = norm(n);

    W = acos ( dot( n, I ) / nnorm );
    if dot(n,J) < 0
        W = 2 * pi - W;
    end

    %% Argument of perigee
    
    % Vector with direction corresponding to apse line
    apse = cross(v,h) / mu - r / rnorm;
    apsenorm = norm(apse);
    
    w = acos ( dot(n, apse)  / nnorm / apsenorm );
    if dot(apse,K) < 0
        w = 2 * pi - w;
    end

    %% True anomaly
    theta = acos ( dot(r,apse) / rnorm / apsenorm );
    if dot(r,v) < 0
        theta = 2 * pi - theta;
    end

    oe_vec(:,ii) = [
        a;
        e;
        rad2deg(i);
        rad2deg(w);
        rad2deg(W);
        rad2deg(theta);
        ];
end
end