function [tspan, oe_vec, ss_vec, geo_vec] = propagator02_ODE_DECHAMPS_FAYT (oe0, tspan, mu)
% This function provides an orbital propagation assuming Keplerian motion
% under the two-body assumption with J2 perturbation.
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
%   - mu        : Gravitational body parameter      [km^3/s^2]
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

% Setting of the solver option
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% Numerical integration through ODE45 solver
[~, ss_vec] = ode45( @(t,ss_vec) keplereq3D(t, ss_vec, mu), ...
    tspan, ss0, options);

% Transformation to orbital elements
oe_vec = cart2kepl_KZ(ss_vec', mu)';

geo_vec = zeros(length(ss_vec),3);

for i = 1:length(ss_vec)
     [h, lambda, phi] = ecef2geodetic([ss_vec(i,1); ss_vec(i,2); ss_vec(i,3)]);
     geo_vec(i,:) = [h rad2deg(lambda) rad2deg(phi)];
end    

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
ddt(1:3) = data(4:end);         % [ x1 y1 z1 ] = [xdot ydot zdot]

% Second derivatives (accelerations) come from force(s) in presence
% Computation of acceleration field
acceleration = accel_field(data(1:3),mu);

% Storage of second derivatives from acceleration field
ddt(4:end) = acceleration(:);

end


function A = accel_field(vec, mu)
% This function provides the acceleration field responsible for the
% movement of the object in the system. The acceleration is given by the
% force in presence. The force field is given by the gradient of the
% potential field. The potential field is given by the spherical harmonics
% expansion truncated to the first term (J2)
% 
%           U = mu / r * (1 - J2 / 2 * R^2 / r^2 * P_2[sin(phi)]
%                               F = grad(U)
% 
% Where R is the average Earth radius, J2 the adimensional perturbation
% coefficient and P_2[sin(phi)] is the second Legendre polynomial : 
%                   P_2[sin(phi)] = 3 * sin(phi)^2 - 1
%
% INPUTS
%   - cart_vec  : Cartesian position vector (1x3)       [m]
%   - mu        : Earth gravitational parameter         [m^3/s^2]
%
% OUTPUTS
%   - A         : Cartesian acceleration field (1x3)    [m/s^2]

% Constants
R = 6371900;           % Earth's average radius (UGGI)     [m]
J = 1.082629e-3;       % Adimensional J2 term              [-]

% Spherical coordinates
r = norm(vec);                      % Radius                [m]
rho = sqrt(vec(1)^2 + vec(2)^2);    % Cylindrical radius    [m]
phi = atan2(rho, vec(3));           % Azimuthal angle       [rad]
the = atan2(vec(2),vec(1));         % Polar angle           [rad]

% Direction cosines
sp = sin(phi); cp = cos(phi);
st = sin(the); ct = cos(the);

% Gradient of potential field could have been computed as :
% syms rv phiv
% U = mu / rv * (1 - J2 / 2 * R^2 / rv^2 * (3 * sin(phiv)^2 - 1)
% accel = [ diff(U,rv); 1/rv * diff(U,phiv); 0];
% A_sph = subs(accel, [rv phiv], [r phi]);
% But the symbolic variable implementation is time-consuming so derivation
% has been made by hand and directly given for computation by the code


% Acceleration in radial direction
% accel_rad = mu  / 2 / r^4 * (   ...
%       9 * J * R^2 * sp^2       ...
%     - 3 * J * R^2              ...
%     - 2 * r^2               );

accel_rad = - 3/2 * J * mu / r^2 * (R/r)^2 * (3 * cp^2 - 1);


% Acceleration in azimuthal direction
% accel_azi =  3 * mu * J * R^2 * sp * cp / r^4;
accel_azi = -3/2 * J * mu / r^2 * (R/r)^2 * sp * cp;

% Acceleration field in spherical coordinates
A_sph = [ 
    accel_rad ; 
    accel_azi ;
            0 ;
         ];


% Transformation matrix from spherical to cartesian
M_sph2cart = [ 
    sp * ct    cp * ct    -st;
    sp * st    cp * st     ct;
         cp        -sp      0;
         ];
     

% Acceleration field in cartesian coordinates
A = M_sph2cart * A_sph;

A = - mu / r^3 * [vec(1); vec(2); vec(3)]  ...
+ 3/2 * J * mu * R^2 / r^4 * [...
    vec(1)/r * (5 * vec(3)^2 / r^2 - 1) ;
    vec(2)/r * (5 * vec(3)^2 / r^2 - 1) ;
    vec(3)/r * (5 * vec(3)^2 / r^2 - 3) ;];

end



function [h, lambda, phi] = ecef2geodetic(rECEF)

% 
% [h, lambda, phi] = ecef2geodetic(rECEF)
% 
% Provide the geodetic coordinates of a position vector in the ECEF frame.
% 
% Inputs:
%   - rECEF = 3-element column vector defining a position in the ECEF frame
%            [m]
% 
% Outputs:
%   - h = altitude from the reference ellipsoid [m]
%   - lambda = geodetic longitude [rad]
%   - phi = geodetic latitude [rad]
% 
% Ref: MATLAB function ecef2geodetic.
% 
% Lamberto Dell'Elce
% 

persistent a e2 ep2 f b

if isempty(a),
    % Ellipsoid constants
    a  = 6378.137e3; % Semimajor axis
    e2 = 0.081819190842^2; % Square of first eccentricity
    ep2 = e2 / (1 - e2); % Square of second eccentricity
    f = 1 - sqrt(1 - e2); % Flattening
    b = a * (1 - f); % Semiminor axis
end

x = rECEF(1);
y = rECEF(2);
z = rECEF(3);

%% Longitude
lambda = atan2(y,x); % [rad]

%% Latitude
% Distance from Z-axis
RHO = hypot(x,y);

% Bowring's formula for initial parametric (beta) and geodetic (phi)
% latitudes
beta = atan2(z, (1 - f) * RHO);
phi = atan2(z   + b * ep2 * sin(beta)^3,...
            RHO - a * e2  * cos(beta)^3);

% Fixed-point iteration with Bowring's formula
% (typically converges within two or three iterations)
betaNew = atan2((1 - f)*sin(phi), cos(phi));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    phi = atan2(z + b * ep2 * sin(beta)^3, ...
        RHO - a * e2  * cos(beta).^3); % [rad]
    betaNew = atan2((1 - f)*sin(phi), cos(phi));
    count = count + 1;
end

%% Altitude
% Calculate ellipsoidal height from the final value for latitude
sin_phi = sin(phi);
N = a / sqrt(1 - e2 * sin_phi^2);
h = RHO * cos(phi) + (z + e2 * N * sin_phi) * sin_phi - N; % [m]

end





