function [tspan, oe_vec, ss_vec] = propagator03_ODE_DECHAMPS_FAYT_drag_only (oe0, tspan, mu, ISS_prop)
% This function provides an orbital propagation assuming Keplerian motion
% under the two-body assumption with the atmospheric drag perturbation.
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
%   - ISS_prop  : Array of the ISS properties       [various]
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
[~, ss_vec] = ode45( @(t,ss_vec) keplereq3D(t, ss_vec, mu, ISS_prop), ...
    tspan, ss0, options);

% Transformation to orbital elements
oe_vec = cart2kepl_KZ(ss_vec', mu)';

geo_vec = zeros(length(ss_vec),3);

for i = 1:length(ss_vec)
     [h, lambda, phi] = ecef2geodetic([ss_vec(i,1); ss_vec(i,2); ss_vec(i,3)]);
     geo_vec(i,:) = [h rad2deg(lambda) rad2deg(phi)];
end    

end

function ddt = keplereq3D(~, data, mu, ISS_prop)
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
acceleration = accel_field(data, mu, ISS_prop);

% Storage of second derivatives from acceleration field
ddt(4:end) = acceleration(:);

end


function A = accel_field(ss_vec, mu, ISS_prop)
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
% The atmospheric drag is modeled by an additional force
%                fA = -.5 * Cd * Aref / m * rho * vTAS^2

%
% INPUTS
%   - cart_vec  : State space vector (1x6)              [m & m/s]
%   - mu        : Earth gravitational parameter         [m^3/s^2]
%
% OUTPUTS
%   - A         : Cartesian acceleration field (1x3)    [m/s^2]

% Constants
R = 6371900;            % Earth's average radius (UGGI)     [m]
% J2 = 1.082629e-3;       % Adimensional J2 term              [-]
J = 1.082629e-3;       % Adimensional J2 term              [-]

% Spherical coordinates
r = norm(ss_vec(1:3));                  % Radius                [m]
rho_r = sqrt(ss_vec(1)^2 + ss_vec(2)^2);  % Cylindrical radius    [m]
phi = atan2(rho_r, ss_vec(3));            % Azimuthal angle       [rad]
the = atan2(ss_vec(2),ss_vec(1));       % Polar angle           [rad]

% Direction cosines
sp = sin(phi); cp = cos(phi);
st = sin(the); ct = cos(the);

% Gradient of potential field could have been computed following :
% syms r_var phi_var
% U = mu / r_var * (1 - J2 / 2 * R^2 / r_var^2 * (3 * sin(phi_var)^2 - 1)
% accel = [ diff(U,r); 1/r * diff(U,phi); 0];
% A_sph = subs(accel, [r_var phi_var], [r phi]);
% But the symbolic variable implementation is time-consuming so derivation
% has been made by hand and directly given for computation by the code


% Acceleration in radial direction
% accel_rad = mu  / 2 / r^4 * (   ...
%       9 * J2 * R^2 * sp^2       ...
%     - 3 * J2 * R^2              ...
%     - 2 * r^2               );

accel_rad = - 3/2 * J * mu / r^2 * (R/r)^2 * (3 * cp^2 - 1);

% Acceleration in azimuthal direction
% accel_azi =  3 * mu * J2 * R^2 * sp * cp / r^4;

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

A = - mu / r^3 * [ss_vec(1); ss_vec(2); ss_vec(3)];  
% + 3/2 * J * mu * R^2 / r^4 * [...
%     ss_vec(1)/r * (5 * ss_vec(3)^2 / r^2 - 1) ;
%     ss_vec(2)/r * (5 * ss_vec(3)^2 / r^2 - 1) ;
%     ss_vec(3)/r * (5 * ss_vec(3)^2 / r^2 - 3) ;];





m_ISS = ISS_prop(1);
Cd_ISS = ISS_prop(2);
S_ref_ISS = ISS_prop(3);

vel_vec = ss_vec(4:end);
rdot = norm(vel_vec);
vel_vec_unitary = vel_vec / rdot;

earth_ang_vel = 2*pi/86400;
earth_ang_vel = [0, 0, earth_ang_vel];

% If the Matlab command "cross" is applied to vectors, they must have a 
% length of 3.

vtas = rdot + cross(ss_vec(1:3), earth_ang_vel); 

rho_atm = harris_priester(ss_vec(1:3));

f_atm_norm = 0.5 * Cd_ISS * S_ref_ISS * rho_atm * vtas.^2 / m_ISS;
f_atm = - f_atm_norm * vel_vec_unitary;

A = A + f_atm*1.75e2;

end




function rho = harris_priester(pos_vec)


% Initialise Harris-Priester data
HP_tab = harris_priester_init;

R = 6371900;            % Earth's average radius (UGGI)     [m]

h = norm(pos_vec) - R;  % Satellite's height            [m]

i = 1;

if h < HP_tab(1,1)
    error('Altitude too low');
else
    while HP_tab(i,1) < h
        i = i + 1;        
    end    
end

hi = HP_tab(i-1,1);
hiplus = HP_tab(i,1);

rho_m_hi = HP_tab(i-1,2);
rho_m_hiplus = HP_tab(i,2);

rho_M_hi = HP_tab(i-1,3);
rho_M_hiplus = HP_tab(i,3);

Hm = (hi - hiplus) / log( rho_m_hiplus / rho_m_hi );
HM = (hi - hiplus) / log( rho_M_hiplus / rho_M_hi );

rho_m_h = rho_m_hi * exp( (hi - h) / Hm);
rho_M_h = rho_M_hi * exp( (hi - h) / HM);

n = 2 ; % Low inclination orbit

e_r = pos_vec / norm(pos_vec);

% delta = deg2rad(11 + 30 / 3600);
delta = deg2rad(23.43929111);
alpha = deg2rad(20);
lambda = deg2rad(30);

e_b = [
    cos(delta) * cos(alpha + lambda);
    cos(delta) * sin(alpha + lambda);
    sin(delta);
    ];

cos_psi = (.5 * ( dot( e_r , e_b ) + 1 ) ) ^ (n/2);
rho = rho_m_h + cos_psi * (rho_M_h - rho_m_h);

end

function [tab] = harris_priester_init()

tab = [...
    100,   497400,  497400;
    120,    24900,   24900;
    130,     8377,    8710;
    140,     3899,    4059;
    150,     2122,    2215;
    160,     1263,    1344;
    170,    800.8,   875.8;
    180,    528.3,     601;   
    190,    361.7,   429.7;
    200,    255.7,   316.2;
    210,    183.9,   239.6;
    220,    134.1,   185.3;
    230,    99.49,   145.5;
    240,    74.88,   115.7;
    250,    57.09,   93.08;
    260,    44.03,   75.55;
    270,    34.30,   61.82;
    280,    26.97,   50.95;
    290,    21.39,   42.26;
    300,    17.08,   35.26;
    320,    10.99,   25.11;
    340,    7.214,   18.19;
    360,    4.824,   13.37;
    380,    3.274,   9.955;
    400,    2.249,   7.492;
    420,    1.558,   5.684;
    440,    1.091,   4.355;
    460,    .7701,   3.362;
    480,    .5474,   2.612;
    500,    .3916,   2.042;
    520,    .2819,   1.605;
    540,    .2042,   1.267;
    560,    .1488,   1.005;
    580,    .1092,   .7997;
    600,    .0807,    .639;
    620,   .06012,   .5123;
    640,   .04519,   .4121;
    660,    .0343,   .3325;
    680,   .02632,   .2691;
    700,   .02043,   .2185;
    720,   .01607,   .1779;
    730,   .01281,   .1452;
    740,   .01036,    .119;   
    780,  .008496,  .09776;
    800,  .007069,  .08059;
    840,   .00468   .05741;
    880,    .0032,   .0421;
    920,   .00221,   .0313;
    960,   .00156,   .0236;
    1000,  .00115,   .0181;
];


tab(:,1) = tab(:,1) * 1e3;
tab(:,2:3) = tab(:,2:3) * 1e-12;

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