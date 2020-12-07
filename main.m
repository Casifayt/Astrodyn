%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AERO0024 - ASTRODYNAMICS %%%
%%%%%% ULiege - 2020-2021 %%%%%%

%%% Authors
% Axel DECHAMPS - S164160
% Casimir  FAYT - S196244


%%%%%% Orbital propagator %%%%%%
% This project consists in the development of an orbital propagator 
% of increasing complexity. The central body is the Earth.

close all; clear; clc; format long;

exo = input(['Please select exercise :\n' ...
    ' 1 for DF1 propagator (two-body model)\n' ...
    ' 2 for DF2 propagator (with J2 perturbation)\n' ...
    ' 3 for DF3 propagator (with J2 and atmospheric drag) \n'...
    ' 4 for comparison with VENuS satellite\n']);

%% Constants
mu = 398600.4418e9;         % Earth gravitational parameter [m^3/s^2]     
ISS_m = 410500;             % ISS's mass                    [kg]
ISS_Cd = 2;                 % ISS's drag coefficient        [-]
ISS_A = 1641;               % ISS's area                    [m^2]

ISS_prop = [ISS_m, ISS_Cd, ISS_A];

MATLABc = {
    [0, 0.4470, 0.7410];
    [0.8500, 0.3250, 0.0980];
    [0.9290, 0.6940, 0.1250];
    [0.4940, 0.1840, 0.5560];
    [0.4660, 0.6740, 0.1880];
    [0.3010, 0.7450, 0.9330];
    [0.6350, 0.0780, 0.1840];
    };


% Time properties
tmax = 86400;
N = 1e3;
dt = tmax/N;
tspan = 0:dt:tmax;

var2clear = {'exo', 'mu', 'ISS_m', 'ISS_Cd', 'ISS_A', 'ISS_prop', 'MATLABc', ...
     'N', 'dt', 'tspan'};
 
%% Initial ISS orbital parameters
e_ISS = .001379;            % Eccentricity      [-]
a_ISS = 6794.57e3;          % Semi-major axis   [m]
i_ISSd = 51.6445;           % Inclination       [deg]
i_ISSr = deg2rad(i_ISSd);   % Inclination       [rad]
W_ISSd = 128.8777;          % RAAN              [deg]
W_ISSr = deg2rad(W_ISSd);   % RAAN              [rad]
w_ISSd = 25.5173;           % Perigee argument  [deg]
w_ISSr = deg2rad(w_ISSd);   % Perigee argument  [rad]
M_ISSd = 146.2321;          % Mean anomaly      [deg]
M_ISSr = deg2rad(M_ISSd);   % Mean anomaly      [rad]


% True anomaly from mean anomaly
% Fourier expansion from wikipedia
% See https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
theta_ISSr = M_ISSr + ...
    (2*e_ISS - e_ISS^3/4) * sin(M_ISSr)   ...
    + 5/4   * e_ISS^2 * sin(2 * M_ISSr)   ...
    + 13/12 * e_ISS^3 * sin(3 * M_ISSr);

theta_ISSd = rad2deg(theta_ISSr);

% Initial orbital elements of ISS
% Degrees
oe_ISSd = [a_ISS, e_ISS, i_ISSd, w_ISSd, W_ISSd, theta_ISSd];
% Radians
oe_ISSr = [a_ISS, e_ISS, i_ISSr, w_ISSr, W_ISSr, theta_ISSr];

% Relative tolerance setting for ode45 solver
% This setting has been optimised by a convergence study (see report for
% more informations, and report_plots/report_plots.m for the script
reltol = 1e-11;

% This cell array allows to clean the workspace by adding the name of
% unuseful variables into it, destined to be wiped out of the memory at the
% end of the code such that the workspace is crystal clear with only
% variables of interest remaining
var2clear = [var2clear, 'a_ISS', 'e_ISS', 'i_ISSd', 'i_ISSr', 'W_ISSd', ...
    'W_ISSr', 'w_ISSd', 'w_ISSr', 'M_ISSd', 'M_ISSr', 'theta_ISSr', ...
    'theta_ISSd', 'oe_ISSd', 'oe_ISSr', 'reltol'];

if exo ~= 4
    fprintf(['\nComputation of the DF' num2str(exo) ...
        ' propagator for a \n '  num2str(tmax/3600) ...
        'h propagation with a time step of ' num2str(dt) 's\n\n']);
end

if exo == 1
    %% Two-body propagator %%
    % S3L orbital propagator

    [~, oe_S3L, ~, ce_S3L] = orbprop(...
        oe_ISSd,                ...     % Initial ISS orbital elements
        'time',     tmax,       ...     % Total time of propagation
        'dt',       dt,         ...     % Time step
        'fmodel',   [0 0 0 0 0] );      % No perturbation

    % DF1 orbital propagator
    [~, oe_DF1, ce_DF1]  =  propagator01_DF(...
        oe_ISSr, ...    % Initial ISS orbital elements
        tspan,   ...    % Vector of time properties
        mu,      ...    % Central body gravitational parameter
        reltol   );     % Relative tolerance for ode45 solver

    % Analytical propagation using Kepler equation
    oe_ISS_KEP = oe_ISSr;
    oe_ISS_KEP(end) = M_ISSr;
    [~, oe_KEP, ce_KEP] = analytical_propagator01(...
        oe_ISS_KEP, ...  % Initial ISS orbital elements with mean anomaly
        tspan,      ...  % Vector of time properties
        mu          );   % Central body gravitational parameter
    
    % Plots comparisons
    vec_cell = { oe_DF1; 'DF1 propagator';
        oe_S3L; 'S3L propagator';
        % Comment next line if analytical propagator not wanted in plot
        oe_KEP; 'Kepler analytical propagator';
        };
    
    coordinates_comparison(vec_cell, tspan, MATLABc, 'Keplerian');

    % Ground tracks
    f = figure;
    f.Name = ('Ground tracks (fixed Earth)');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_DF1, 'DF1 propagator', 0, dt);
    subplot(2,1,2);
    grdtrk(ce_S3L, 'S3L propagator', 0, dt);
    
    f = figure;
    f.Name = ('Ground tracks (rotating Earth)');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_DF1, 'DF1 propagator', 1, dt);
    subplot(2,1,2);
    grdtrk(ce_S3L, 'S3L propagator', 1, dt);
    

    orb_trk_3d_DF(ce_DF1(:,1:3));
    
    % Comparison of final elements
    vec_cell = {ce_DF1, oe_DF1, 'DF1', ...
                ce_S3L, oe_S3L, 'S3L'};
    final_elements_print(vec_cell);
    
    fprintf('\n\n\n');
    
    vec_cell = {ce_KEP, oe_KEP, 'KEPL', ...
                ce_DF1, oe_DF1, 'DF1'};
    
    final_elements_print(vec_cell);
    
    var2clear = [var2clear, 'oe_ISS_KEP', 'f', 'vec_cell' ];
    
elseif exo == 2

%% J2-term (Earth's oblateness) %%

    % S3L orbital propagator
    [~, oe_S3L, ~, ce_S3L] = orbprop(...
        oe_ISSd,                 ...    % Initial ISS orbital elements
        'time',     tspan(end),  ...    % Total time of propagation
        'dt',       dt,          ...    % Time step
        'fmodel',   [1 0 0 0 0]  );     % J2 perturbation
        

    % DF2 orbital propagator
    [~, oe_DF2, ce_DF2]  =  propagator02_DF(...
        oe_ISSr, ...    % Initial ISS orbital elements
        tspan,   ...    % Vector of time properties
        mu,      ...    % Central body gravitational parameter
        reltol   );     % Relative tolerance setting for ode45 solver

    
    % Analytical propagation using Gauss perturbation equations
    oe_ISS_GAUSS = oe_ISSr;
    oe_ISS_GAUSS(end) = M_ISSr;
    [~, oe_GAUSS, ce_GAUSS] = analytical_propagator02(...
        oe_ISS_GAUSS, ... % Initial ISS orbital elements with mean anomaly
        tspan,        ... % Vector of time properties
        mu            );  % Central body gravitational parameter

    
    % Plots comparisons
    vec_cell = {
        oe_DF2; 'DF2 propagator';
        oe_S3L; 'S3L propagator';
        % Comment next line if analytical propagator not wanted in plot
        oe_GAUSS; 'Gauss analytical propagator';
        };
    
    coordinates_comparison(vec_cell, tspan, MATLABc, 'Keplerian');
    
    % Ground tracks
    f = figure;
    f.Name = ('Ground tracks of the orbit');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_DF2(:,1:3), 'Fixed Earth', 0, dt);
    subplot(2,1,2);
    grdtrk(ce_DF2(:,1:3), 'Rotating Earth', 1, dt);
    
    orb_trk_3d_DF(ce_DF2(:,1:3));
    
    % Comparison of final elements
    vec_cell = {ce_DF2, oe_DF2, 'DF2', ...
                ce_S3L, oe_S3L, 'S3L', ...
        };
    
    final_elements_print(vec_cell);
    
    
    var2clear = [var2clear, 'oe_ISS_GAUSS', 'f', 'vec_cell'];

    
elseif exo == 3
%% Earth's Atmosphere %%

    % S3L orbital propagator
    [~, oe_S3L, ~, ce_S3L] = orbprop(...
        oe_ISSd,                ...     % Initial ISS orbital elements
        'time',     tmax,       ...     % Total time of propagation
        'dt',       dt,         ...     % Time step
        'fmodel',   [1 1 0 0 0] );      % J2 and drag perturbations

    % S3Ld orbital propagator (S3L with drag only)
    [~, oe_S3Ld, ~, ce_S3Ld] = orbprop(...
        oe_ISSd,                ...     % Initial ISS orbital elements
        'time',   tmax,         ...     % Total time of propagation
        'dt',     dt,           ...     % Time step
        'fmodel', [0 1 0 0 0]   );      % Only drag perturbation    
    
    % DF3 orbital propagator
    [~, oe_DF3, ce_DF3]  =  propagator03_DF(oe_ISSr, tspan, mu, ISS_prop);
    
    % DF3d orbital propagator (DF3 with drag only)
    [~, oe_DF3d, ce_DF3d]  =  propagator03_DF_drag_only(...
        oe_ISSr,  ...   % Initial ISS orbital elements
        tspan,    ...   % Vector of time properties
        mu,       ...   % Central body gravitational parameter
        ISS_prop  );    % Properties of ISS satellite
    
    [a_red_anal,a_red_num,a_red_S3L] = analytical_propagator03(...
        oe_DF3d(:,1),  ...  % Semi-major axis evolution by DF3d
        oe_S3Ld(:,1),  ...  % Semi-major axis evolution by S3Ld
        tspan       );
    % Plots comparisons
    vec_cell = {
        oe_DF3; 'DF3 propagator';
        oe_S3L; 'S3L propagator';
        };
    
    coordinates_comparison(vec_cell, tspan, MATLABc, 'Keplerian');

    % Ground tracks
    f = figure;
    f.Name = ('Ground tracks');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_DF3, 'Fixed Earth',0,dt);
    subplot(2,1,2);
    grdtrk(ce_S3L, 'Rotating Earth',1,dt);
    
    orb_trk_3d_DF(ce_DF3(:,1:3));
    
    % Final elements print
    vec_cell = {ce_DF3, oe_DF3, 'DF3', ...
                ce_S3L, oe_S3L, 'S3L', ...
    };

    final_elements_print(vec_cell);
    
    vec_cell = {ce_DF3d, oe_DF3d, 'DF3d', ...
                ce_S3L, oe_S3L, 'S3L', ...
    };

    final_elements_print(vec_cell);

    var2clear = [var2clear, 'a_red_anal', 'a_red_num', 'a_red_S3L', ...
        'f', 'vec_cell'];

elseif exo == 4
%% Comparison with actual satellite data : Ven탎 %%
    
    fprintf(['\nComparison of the propagators using the \n' ...
        'VEN킪 satellite data\n\n']);

    % Time properties
    tmax = 89175;                                   % Time between the two dates of the two TLEs [s]
    N = 1e3;
    dt = tmax/N;
    tspan = 0:dt:tmax;
    
    % Conversion of the TLE into the ECI frame
    [oe_start, oe_final, ss_start, ss_final] = TLE2ECI_VENUS(mu);
    
    % Initial Ven탎 Orbital Parameters
    a_VNS = oe_start(1);                    % Semi-major axis   [m]
    e_VNS = oe_start(2);                    % Eccentricity      [-]
    i_VNSd = oe_start(3);                   % Inclination       [deg]
    i_VNSr = deg2rad(i_VNSd);               % Inclination       [rad]
    w_VNSd = oe_start(4);                   % Perigee argument  [deg]
    w_VNSr = deg2rad(w_VNSd);               % Perigee argument  [rad]
    W_VNSd = oe_start(5);                   % RAAN              [deg]
    W_VNSr =  deg2rad(W_VNSd);              % RAAN              [rad]
    theta_VNSd = oe_start(6);               % True anomaly      [deg]
    theta_VNSr = deg2rad(theta_VNSd);       % True anomaly      [rad]

    oe_VNSd = [a_VNS, e_VNS, i_VNSd, w_VNSd, W_VNSd, theta_VNSd];
    oe_VNSr = [a_VNS, e_VNS, i_VNSr, w_VNSr, W_VNSr, theta_VNSr];

    % Ven탎 ballistic properties

    m_VNS = 265;                            % Mass of Ven탎 [kg]
    Cd_VNS = 2;                             % Drag coefficient [-]
    A_VNS = 2.2619;                         % Cross-section of Venus [m^2]

    VNS_prop = [m_VNS, Cd_VNS, A_VNS];

    
    % S3L orbital propagator
    [~, oe_S3L, ~, ce_S3L] = orbprop(...
        oe_VNSd,               ...  % Initial Ven탎 orbital elements
        'time',   tmax,        ...  % Total time of propagation
        'dt',     dt,          ...  % Time step
        'fmodel', [1 1 0 0 0]  );   % J2 and drag perturbations

    % Propagation of Ven탎 using the DF3 propagator
    [~, oe_DF3, ce_DF3]  =  propagator03_DF(...
        oe_VNSr,  ...       % Initial Ven탎 orbital elements
        tspan,    ...       % Vector of time properties
        mu,       ...       % Central body gravitational parameter
        VNS_prop  );        % Properties of Ven탎 satellite
    
    % Propagation of Ven탎 using the SGP4 propagator
    [oe_SGP4,ce_SGP4] = SGP4_VENUS(mu, tspan);
    
    % Plots comparisons
    vec_cell = {
        oe_DF3;  'DF3 propagator';
        oe_S3L;  'S3L propagator';
        oe_SGP4; 'SGP4 propagator';
        };
    
    coordinates_comparison(vec_cell, tspan, MATLABc, 'Keplerian');

    % Ground track
    f = figure;
    f.Name = 'Ground tracks';
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_DF3, 'Fixed Earth',0,dt);
    subplot(2,1,2);
    grdtrk(ce_DF3, 'Rotating Earth',1,dt);
    
    orb_trk_3d_DF(ce_DF3(:,1:3));
    
    vec_cell = {ce_DF3,   oe_DF3,   'DF3',  ...
                ce_SGP4,  oe_SGP4,  'SGP4' };
    final_elements_print(vec_cell);
    
    vec_cell = {ce_DF3,   oe_DF3,   'DF3',  ...
                ss_final', oe_final', 'CLS'   };
    final_elements_print(vec_cell);
    
        var2clear = [var2clear, 'dt', 'tspan', 'f', 'vec_cell', ...
            'a_VNS', 'e_VNS', 'i_VNSd', 'i_VNSr', 'w_VNSr', 'w_VNSd', ...
            'W_VNSd', 'W_VNSr', 'theta_VNSr', 'theta_VNSd', 'oe_VNSd', 'oe_VNSr', ...
            'VNS_prop', 'm_VNS', 'A_VNS', 'Cd_VNS', 'oe_final', 'oe_start',...
            'ss_final', 'ss_start'];
    
end

% Unuseful variables clearing
var2clear = [var2clear, 'var2clear'];
clear(var2clear{:});



%% Other functions
function final_elements_print(vec_cell)
% This function prints final coordinates of up to 3 vectors and computes
% the relative least-squared error about the coordinates.
% INPUTS
%   - vec_cell : Cell of the vectors to print and their title, ordered as
%       - vec_cell{3i+1} - Cartesian coordinates of first vector
%       - vec_cell{3i+2} - Keplerian coordinates of first vector
%       - vec_cell{3i+3} - Title fo first vector


    vec1_ce = vec_cell{1}; vec1_oe = vec_cell{2}; title1 = vec_cell{3};    
    vec2_ce = vec_cell{4}; vec2_oe = vec_cell{5}; title2 = vec_cell{6};

	fprintf([ '\n --------------------------------------------' ...
        '\n| COORDINATES COMPARISON BETWEEN ' title1 ' AND ' title2 '|' ...
        '\n --------------------------------------------']);
    
    % Cartesian print
    x2 = vec2_ce(end,1); xdot2 = vec2_ce(end,4);
    y2 = vec2_ce(end,2); ydot2 = vec2_ce(end,5);
    z2 = vec2_ce(end,3); zdot2 = vec2_ce(end,6);
    
    x1 = vec1_ce(end,1); xdot1 = vec1_ce(end,4);
    y1 = vec1_ce(end,2); ydot1 = vec1_ce(end,5);
    z1 = vec1_ce(end,3); zdot1 = vec1_ce(end,6);
    
    erel_x = abs( (x2 - x1) / x2);
    erel_y = abs( (y2 - y1) / y2);
    erel_z = abs( (z2 - z1) / z2);
    
    erel_xdot = abs( (xdot2 - xdot1) / xdot2);
    erel_ydot = abs( (ydot2 - ydot1) / ydot2);
    erel_zdot = abs( (zdot2 - zdot1) / zdot2);
    
    fprintf(['\n1. Final cartesian coordinates\n'             ...
       '                  ' title1 '           ' title2 '      Error\n'       ...
       'x [km] =         %.2f     %.2f    %.1e %%\n'           ...
       'y [km] =         %.2f      %.2f     %.1e %%\n'         ...
       'z [km] =         %.2f     %.2f      %.1e %%\n'         ...
       'xdot [km/s] =    %.2f         %.2f      %.1e %%\n'     ...
       'ydot [km/s] =    %.2f         %.2f      %.1e %%\n'     ...
       'zdot [km/s] =    %.2f           %.2f      %.1e %%\n'], ...
        x1/1000,    x2/1000, 100 *    erel_x,            ...
        y1/1000,    y2/1000, 100 *    erel_y,            ...
        z1/1000,    z2/1000, 100 *    erel_z,            ...
     xdot1/1000, xdot2/1000, 100 * erel_xdot,            ...
     ydot1/1000, ydot2/1000, 100 * erel_ydot,            ...
     zdot1/1000, zdot2/1000, 100 * erel_zdot             ...
     );
    
     fprintf(['Least squared error on cartesian vectors : \n' ...
         'On position : %.2e %%\n'                            ...
         'On velocity : %.2e %%\n']                         , ...
         sqrt(erel_x^2 + erel_y^2 + erel_z^2)              , ...
         sqrt(erel_xdot^2 + erel_ydot^2 + erel_zdot^2)       ...
     );
 
    % Keplerian print
    a2 = vec2_oe(end,1); w2 = vec2_oe(end,4);
    e2 = vec2_oe(end,2); W2 = vec2_oe(end,5);
    i2 = vec2_oe(end,3); t2 = vec2_oe(end,6);
    
    a1 = vec1_oe(end,1); w1 = vec1_oe(end,4);
    e1 = vec1_oe(end,2); W1 = vec1_oe(end,5);
    i1 = vec1_oe(end,3); t1 = vec1_oe(end,6);
    
    erel_a = abs(a2 - a1) / a2;
    erel_e = abs(e2 - e1) / e2;
    erel_i = abs(i2 - i1) / i2;
    erel_w = abs(w2 - w1) / w2;
    erel_W = abs(W2 - W1) / W2;
    erel_t = abs(t2 - t1) / t2;
    
    fprintf(['\n2. Final Keplerian coordinates \n'         ...
       '                  ' title1 '         ' title2 '     Error\n'       ...
       'a [km] =        %.2f    %.2f    %.1e %%\n'          ...
       'e [-] =         %.2e  %.2e    %.1e %%\n'            ...
       'i [deg] =       %.2f      %.2f      %.1e %%\n'      ...
       '\x03C9 [deg] =       %.2f      %.2f      %.1e %%\n' ...
       '\x03A9 [deg] =       %.2f     %.2f     %.1e %%\n'   ...
       '\x03B8 [deg] =       %.2f     %.2f     %.1e %%\n'], ...
     a1/1000, a2/1000, 100 * erel_a,                  ...
          e1,      e2, 100 * erel_e,                  ...
          i1,      i2, 100 * erel_i,                  ...
          w1,      w2, 100 * erel_w,                  ...
          W1,      W2, 100 * erel_W,                  ...
          t1,      t2, 100 * erel_t                   ...
     );
end

function coordinates_comparison(vec_cell, tspan, MATLABc, type)
% This function plots a graphical representation of the evolution of the
% elements of a Nx6 coordinates vector across N time steps, given in tspan.
% According to the input type given, it will plot cartesian or keplerian
% coordinates.
% INPUTS
%   - vec_cell : Cell of up to 3 vectors and their title to plot, ordered as
%       - vec_cell{2i+1} - Nx6 vector of elements to plot
%       - vec_cell{2i+2} - title of the vector          
%   - tspan     : Vector of time properties             [Nx1]
%   - MATLABc   : List of RGB colors inputs for plots   [Nx3]
%   - type      : 'Keplerian' or 'Cartesian'            [str]

    vec_1 = vec_cell{1};
    title1 = vec_cell{2};

    vec_2 = vec_cell{3};
    title2 = vec_cell{4};

    if numel(vec_cell) > 4
        vec_3 = vec_cell{5};
        title3 = vec_cell{6};
    end

    f = figure;
    
    if strcmp(type, 'Cartesian') == 1 || ...
       strcmp(type, 'cartesian') == 1 || ...
       strcmp(type, 'cart') == 1
        
            f.Name = ('Comparison of Cartesian elements');
            f.WindowState = 'maximized';

            subplot(3,2,1)
            plot( tspan/3600 ,  vec_1(:,1)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,1)/1000 , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,1)/1000 , 'Color' , MATLABc{3}); 
            end
            title('x'); ylabel('Position [km]'); xlabel('Time [hours]');

            subplot(3,2,3)
            plot( tspan/3600 ,  vec_1(:,2)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,2)/1000 , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,2)/1000 , 'Color' , MATLABc{3}); 
            end
            title('y'); ylabel('Position [km]'); xlabel('Time [hours]');

            subplot(3,2,5)
            plot( tspan/3600 ,  vec_1(:,3)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,3)/1000 , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,3)/1000 , 'Color' , MATLABc{3}); 
            end
            title('z'); ylabel('Position [km]'); xlabel('Time [hours]');

            subplot(3,2,2);
            plot( tspan/3600 ,  vec_1(:,4)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,4)/1000 , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,4)/1000 , 'Color' , MATLABc{3}); 
            end
            title('$\mathbf{\dot{x}}$', 'interpreter', 'latex'); 
            ylabel('Velocity [km/s]'); xlabel('Time [hours]');

            subplot(3,2,4);
            plot( tspan/3600 ,  vec_1(:,5)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,5)/1000 , 'Color' , MATLABc{2});
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,5)/1000 , 'Color' , MATLABc{3}); 
            end
            title('$\mathbf{\dot{y}}$', 'interpreter', 'latex'); 
            ylabel('Velocity [km/s]'); xlabel('Time [hours]');

            subplot(3,2,6);
            plot( tspan/3600 ,  vec_1(:,6)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,6)/1000 , 'Color' , MATLABc{2});
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,6)/1000 , 'Color' , MATLABc{3}); 
            end
            title('$\mathbf{\dot{z}}$', 'interpreter', 'latex'); 
            ylabel('Velocity [km/s]'); xlabel('Time [hours]');
            
            if numel(vec_cell) > 4
                leg = legend(title1, title2, title3);
            else
                leg = legend(title1, title2);
            end
            
            set(leg,'Position', [.455 .0125 .125 .075], 'Units', 'normalized');
            leg.Position(1) = 0.5 - leg.Position(3)/2; 
   
    elseif strcmp(type, 'Keplerian') == 1 || ...
       strcmp(type, 'keplerian') == 1 || ...
       strcmp(type, 'kepl') == 1 || strcmp(type, 'Orbital') == 1 || ...
       strcmp(type, 'orbital') == 1 || ...
       strcmp(type, 'orb') == 1
        
            f.Name = ('Comparison of orbital elements');
            f.WindowState = 'maximized';

            subplot(3,2,1);
            plot( tspan/3600 ,  vec_1(:,1)/1000 , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,1)/1000 , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,1)/1000 , 'Color' , MATLABc{3}); 
            end
            title('Semi-major axis'); ylabel('a [km]'); xlabel('Time [hours]');

            subplot(3,2,3);
            plot( tspan/3600 ,  vec_1(:,2) , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,2) , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 ,  vec_3(:,2) , 'Color' , MATLABc{3}); 
            end
            title('Eccentricity'); ylabel('e [-]'); xlabel('Time [hours]');

            subplot(3,2,5);
            plot( tspan/3600 ,  vec_1(:,3) , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,3) , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 , vec_3(:,3) , 'Color' , MATLABc{3}); 
            end
            title('Inclination'); ylabel('i [deg]'); xlabel('Time [hours]');

            subplot(3,2,2);

            for i = 1:length(tspan)
                if vec_1(i,4) > 350
                   vec_1(i,4) = vec_1(i,4) - 360;
                end
                if vec_2(i,4) > 350
                   vec_2(i,4) = vec_2(i,4) - 360;
                end
                if numel(vec_cell) > 4
                    if vec_3(i,4) > 350
                       vec_3(i,4) = vec_3(i,4) - 360;
                    end
                end
            end

            plot( tspan/3600 ,  vec_1(:,4) , 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,4) , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 , vec_3(:,4) , 'Color' , MATLABc{3}); 
            end
            title('Argument of perigee'); 
            ylabel('\omega [deg]'); xlabel('Time [hours]');

            subplot(3,2,4);
            plot( tspan/3600 ,  vec_1(:,5), 'Color' , MATLABc{1}); hold on;
            plot( tspan/3600 ,  vec_2(:,5) , 'Color' , MATLABc{2}); 
            if numel(vec_cell) > 4
               plot( tspan/3600 , vec_3(:,5) , 'Color' , MATLABc{3}); 
            end
            title('RAAN'); ylabel('\Omega [deg]'); xlabel('Time [hours]');

            subplot(3,2,6);
            plot( tspan/3600 ,  vec_1(:,6) , 'Color' , MATLABc{1});  hold on;
            plot( tspan/3600 ,  vec_2(:,6) , 'Color' , MATLABc{2});
            if numel(vec_cell) > 4
               plot( tspan/3600 , vec_3(:,6) , 'Color' , MATLABc{3});
            end
            title('True anomaly'); ylabel('\theta [deg]'); xlabel('Time [hours]');

            if numel(vec_cell) > 4
                leg = legend(title1, title2, title3);
            else
                leg = legend(title1, title2);
            end
            
            set(leg,'Position', [.455 .0125 .125 .075], 'Units', 'normalized');
            leg.Position(1) = 0.5 - leg.Position(3)/2; 
            
    end
end


