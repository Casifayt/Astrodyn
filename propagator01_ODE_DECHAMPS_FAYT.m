function [tspan, oe_vec, ss_vec] = propagator01_ODE_DECHAMPS_FAYT (oe0, tspan, mu, reltol)
% This function provides an orbital propagation assuming Keplerian motion
% under the two-body assumption.
% The EoM are integrated using the ODE45 solver.
% 
% INPUTS
%   - oe0   : Initial vector of keplerian coordinates ordered as :
%       - oe0(1) = a     - semi-major axis          [m]
%       - oe0(2) = e     - orbit eccentricity       [-]
%       - oe0(3) = i     - inclination              [rad]
%       - oe0(4) = omega - argument of perigee      [rad]
%       - oe0(5) = Omega - RAAN                     [rad]
%       - oe0(6) = theta - true anomaly             [rad]
%   - tspan     : Vector of time properties         [s]
%   - mu        : Gravitational body parameter      [m^3/s^2]
%   - reltol    : Relative tolerance of the solver  [-]
% 
% OUTPUTS
%   - tspan     : Vector of time properties         [s]
%   - oe_vec    : Final keplerian coordinates vector (1x6) ordered as the
%   input vector
%   - ss_vec    : Final cartesian coordinates vector (1x6) ordered as :
%       - ss_vec(1:3) = r_vec = [   x    y    z  ]      [m]
%       - ss_vec(4:6) = v_vec = [ xdot ydot zdot ]      [m/s]


% For computational purposes, Cartesian coordinates are used.
ss0 = kepl2cart_KZ(oe0, mu);

% Setting the solver options
options = odeset('RelTol',reltol,'AbsTol',1e-13);

% Numerical integration through ODE45 solver
[~, ss_vec] = ode45( @(t,ss_vec) keplereq3D(t, ss_vec, mu), ...
    tspan, ss0, options);

% Transformation to orbital elements
oe_vec = cart2kepl_KZ(ss_vec', mu)';

end

function ddt = keplereq3D(~, data, mu)
% This function transforms the second-order equation in a system of
% first-order differential equations.
% 
% INPUTS
%   - tspan     : Incremental time step vector      [s]
%   - data      : State-space vector (1x6) ordered as follows
%       - data(1:3) = r_vec = [   x    y    z  ]    [m]
%       - data(4:6) = v_vec = [ xdot ydot zdot ]    [m/s]
%   - mu        : Earth gravitational parameter     [m^3/s^2]
% 
% OUTPUTS
%   - ddt       : Derivatives vector (1x6) ordered as follows
%       - ddt(1:3) = v = [  xdot   ydot   zdot ]    [m/s]
%       - ddt(4:6) = a = [ xddot  yddot  zddot ]    [m/s^2]


% Initialisation of derivatives vector
ddt = zeros(6,1);

% First derivatives (velocities) are given in state-space vector
% Allocation expanded for educational purpose
ddt(1) = data(4);                   % x1 = xdot
ddt(2) = data(5);                   % y1 = ydot
ddt(3) = data(6);                   % z1 = zdot

% Second derivatives (accelerations) come from force(s) in presence
% Computation of acceleration field
acceleration = accel_field(data(1:3), mu);


% Storage of second derivatives from acceleration field
% Allocation expanded for educational purpose
ddt(4) = acceleration(1);           % x2 = xddot
ddt(5) = acceleration(2);           % y2 = yddot
ddt(6) = acceleration(3);           % z2 = zddot

end


function A = accel_field(cart_vec, mu)
% This function provides the acceleration field responsible for the
% movement of the object in the system. The acceleration is given by the
% force in presence. The force field is given by the gradient of the
% potential field. In the two-body approximation, the potential field is
% the gravitational field of the Earth 
% 
%                           U = mu / r
%                          F =  grad(U) 
%
% INPUTS
%   - cart_vec  : Cartesian position vector (1x3)       [m]
%   - mu        : Earth gravitational parameter         [m^3/s^2]
%
% OUTPUTS
%   - A         : Cartesian acceleration field (1x3)    [m/s^2]

% Norm of radius vector
r = norm(cart_vec);

% Acceleration field
A = - mu / r^3 * cart_vec;

end

















