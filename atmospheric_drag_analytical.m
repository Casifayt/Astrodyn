function [a_reduction_anal,a_reduction_num,a_reduction_SL3] = atmospheric_drag_analytical(a_drag,a_SL3) 

% This function computes the reduction of the semi-major axis of the orbit
% assuming Keplerian motion under the two-body assumption with the 
% atmospheric drag perturbation only.

% INPUTS
%   - a_drag                : Semi-major axis computed with the propagator03 [m]
%   - a_SL3                 : Semi-major axis computed with the S3L propagator [m]

% OUTPUTS
%   - a_reduction_anal      : Reduction of the semi-major axis of the orbit [m] (analytical)
%   - a_reduction_num       : Reduction of the semi-major axis of the orbit [m] (numerical - propagator03)
%   - a_reduction_SL3       : Reduction of the semi-major axis of the orbit [m] (numerical - S3L propagator)

%% Parameters

a_i = 6794570;                  % Initial semi-major axis [m]
mu = 398600.4418e9;             % Earth gravitational parameter [m^3/s^2]    
m = 410500;                     % Mass of the ISS [kg]
Cd = 2;                         % Drag coefficient [-]
A = 1641;                       % Cross-section of the ISS [m^2]
t_f = 86400;                    % Final time (24h) [s]
t_i = 0;                        % Initial time [s]
alti = 400;                     % Altitude of the orbit [km]

%% Computation of the Density

rho = harris_priester(alti);    % Atmospheric density [kg/m^3]

%% Computation of the Final Semi-Major Axis

a_f = sqrt(a_i) - (sqrt(mu)*rho*Cd*A)/(2*m)*(t_f-t_i);  % af = sqrt(af)
a_f = a_f^2;                                            % Final semi-major axis [m]

%% Computation of the Reduction in the Semi-Major Axis (Analytical)

diff = a_f - a_i;             % Difference between the initial and final semi-major axes [m]

a_reduction_anal = diff;      % Difference between the initial and final semi-major axes [m] (OUTPUT)

%% Computation of the Reduction in the Semi-Major Axis (Propagator03)

[~,pos] = min(a_drag);

diff = a_drag(pos,1) - a_i;       % Difference between the initial and final semi-major axes [m]

a_reduction_num = diff;           % Difference between the initial and final semi-major axes [m] (OUTPUT)

%% Computation of the Reduction in the Semi-Major Axis (S3L Propagator)

[~,pos] = min(a_SL3);

diff = a_SL3(pos,1) - a_i;        % Difference between the initial and final semi-major axes [m]

a_reduction_SL3 = diff;           % Difference between the initial and final semi-major axes [m] (OUTPUT)

%% Plot

plot([t_i,t_f/3600],[a_i,a_f],'color',[252 186 3]/255,'LineWidth',1.5)
hold on
plot([t_i,t_f/3600],[a_i,a_drag(pos,1)],'color',[15 43 184]/255,'LineWidth',1.5)
hold on
plot([t_i,t_f/3600],[a_i,a_SL3(pos,1)],'color',[201 4 4]/255,'LineWidth',1.5)
h=legend('Analytical Result','Propagator03 (drag only)','S3L Propagator (drag only)');
set(h,'interpreter','Latex','FontSize',12,'Location','south west');
grid minor
xlabel('Time [h]','Fontsize',12)
ylabel('Semi-Major Axis [m]','Fontsize',12)

set(gca,'fontsize',12)
set(gcf, 'position', [300, 200, 700, 500])


end

%% Other Required Functions

function rho = harris_priester(alti)


% Initialise Harris-Priester data
HP_tab = harris_priester_init;

% R = 6371900;            % Earth's average radius (UGGI)     [m]

% h = norm(alti) - R;     % Satellite's height            [m]
h = alti*1e3;   % Satellite's height            [m]

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

e_r = alti / norm(alti);
e_r = [e_r ; e_r ; e_r];

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
