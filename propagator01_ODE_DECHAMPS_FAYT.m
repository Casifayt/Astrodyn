%%%%%% Two-body propagator %%%%%%
% This function provides an orbital propagation assuming Keplerian motion
% The EoM are integrated using the ODE45 solver.

% Governing EoM is
%       ..
%  \vec{ r} = - (mu / r^3) \vec{r}


% INPUTS
%   - vec0      : Initial vector of coordinates

% If Keplerian coordinates are used, elements are ordered as follows : 
%       - vec0(1) = semi-major axis      [m]
%       - vec0(2) = orbit eccentricity   [-]
%       - vec0(3) = inclination          [deg]
%       - vec0(4) = argument of perigee  [deg]
%       - vec0(5) = RAAN                 [deg]
%       - vec0(6) = true anomaly         [deg]

% If Cartesian coordinates are used, elements are ordered as follows : 
%       - vec0(1:3) = Cartesian position [m]
%       - vec0(4:6) = Cartesian velocity [m/s]

%   - tspan     : Time vector properties            [s]
%   - mu        : Gravitational body parameter      [km^3/s^2]

% OUTPUTS
%   - ss_vec    : State vector


function [ss_vec] = propagator01_ODE_DECHAMPS_FAYT (vec0, tspan, mu)

% In the two-body approximation, the orbit is circular.
% All keplerian coordinates are constant except for the true anomaly.
% The evolution of the true anomaly is directly integrated from Kepler's
% laws of motion


% But for computational purposes, Cartesian coordinates are used.
[~, ss_vec] = ode45( @(t, ss_vec) odefcn3D(t ,ss_vec, mu), tspan, vec0);
ss_vec = ss_vec';

end

function drdt = odefcn3D(t, ss_vec, mu)
% This function represents the 2nd order ODE as a system of 1st order 
% ODEs, in three dimensions.

rnorm = norm(ss_vec(1:3));

drdt = [
    odefcn1D(t, [ss_vec(1) ss_vec(4)], mu, rnorm) ;
    odefcn1D(t, [ss_vec(2) ss_vec(5)], mu, rnorm) ;
    odefcn1D(t, [ss_vec(3) ss_vec(6)], mu, rnorm) ;
    ];

end


function ddt = odefcn1D(~, pos, mu, rnorm)
% This function represents the 2nd order ODE as a system of 1st order 
% ODEs, in one dimension.
ddt = zeros(2,1);
ddt(1) = pos(2);
ddt(2) = - mu / rnorm^3 * pos(1);
end


















