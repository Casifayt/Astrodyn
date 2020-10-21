%%%%%% Two-body propagator %%%%%%
% This function provides an orbital propagation assuming Keplerian motion
% The EoM are integrated using the ODE45 solver.
% 
% Governing EoM is
%       ..
%  \vec{ r} = - (mu / r^3) \vec{r}
% 
% 
% INPUTS
%   - oe0   : Initial vector of keplerian coordinates ordered as :
%       - oe0(1) = a -     semi-major axis          [m]
%       - oe0(2) = e -     orbit eccentricity       [-]
%       - oe0(3) = i -     inclination              [rad]
%       - oe0(4) = omega - argument of perigee      [rad]
%       - oe0(5) = Omega - RAAN                     [rad]
%       - oe0(6) = theta - true anomaly             [rad]
%   - tspan     : Vector of time properties         [s]
%   - mu        : Gravitational body parameter      [km^3/s^2]
% 
% OUTPUTS
%   - tspan     : Vector of time properties         [s]
%   - oe_vec    : Final keplerian coordinates vector (1x6) ordered as the
%   input vector
%   - ss_vec    : Final cartesian coordinates vector (1x6) ordered as :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]


function [tspan, oe_vec, ss_vec] = propagator01_ODE_DECHAMPS_FAYT (oe0, tspan, mu)

% In the two-body approximation, the orbit is circular.
% All keplerian coordinates are constant except for the true anomaly.
% The evolution of the true anomaly is directly integrated from Kepler's
% laws of motion

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% But for computational purposes, Cartesian coordinates are used.
ss0 = kepl2cart_KZ(oe0, mu);

[~, ss_vec] = ode45( @(t,ss_vec) keplereq3D(t, ss_vec, mu), ...
    tspan, ss0, options);

oe_vec = cart2kepl_KZ(ss_vec', mu);
oe_vec = oe_vec';
end

function ddt = keplereq3D(~, data, mu)
% This function represents the 2nd order kepler equation of relative motion
% as a system of 1st order ODE, in one dimension.

rnorm = norm(data(1:3));

ddt = zeros(6,1);
ddt(1) = data(4);                   % x1
ddt(2) = data(5);                   % y1
ddt(3) = data(6);                   % z1

ddt(4) = - mu / rnorm^3 * data(1);  % x2
ddt(5) = - mu / rnorm^3 * data(2);  % y2
ddt(6) = - mu / rnorm^3 * data(3);  % z2

end




function ddt = keplereq1D(~, data, mu, rnorm)
% This function represents the 2nd order kepler equation of relative motion
% as a system of 1st order ODE, in one dimension.
ddt = zeros(2,1);
ddt(1) = data(2);
ddt(2) = - mu / rnorm^3 * data(1);
end


















