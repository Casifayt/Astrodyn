%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AERO0024 - ASTRODYNAMICS %%%
%%%%%% ULiège - 2020-2021 %%%%%%

%%% Authors
% Axel DECHAMPS - S164160
% Casimir  FAYT - S196244


%%%%%% Orbital propagator %%%%%%
% This project consists in the development of an orbital propagator 
% of increasing complexity. The central body is the Earth.

close all; clear; clc; format long;

exo = input(['Please select exercise :\n' ...
    ' 1 for two-body\n' ...
    ' 2 for J2 \n' ...
    ' 3 for atmospheric drag (TBD) \n'...
    ' 4 for genuine comparison (TBD) \n']);

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

%% Initial ISS orbital parameters
e_ISS = .001379;                    % Eccentricity      [-]
a_ISS = 6794.57e3;                  % Semi-major axis   [m]
i_ISSd = 51.6445;                   % Inclination       [deg]
i_ISSr = deg2rad(i_ISSd);           % Inclination       [rad]
RAAN_ISSd = 128.8777;               % RAAN              [deg]
RAAN_ISSr =  deg2rad(RAAN_ISSd);    % RAAN              [rad]
omega_ISSd = 25.5173;               % Perigee argument  [deg]
omega_ISSr = deg2rad(omega_ISSd);   % Perigee argument  [rad]
M_ISSd = 146.2321;                  % Mean anomaly      [deg]
M_ISSr = deg2rad(M_ISSd);           % Mean anomaly      [rad]


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
oe_ISSd = [a_ISS, e_ISS, i_ISSd, omega_ISSd, RAAN_ISSd, theta_ISSd];
% Radians
oe_ISSr = [a_ISS, e_ISS, i_ISSr, omega_ISSr, RAAN_ISSr, theta_ISSr];

reltol = 1e-11;

if exo == 1
    %% Two-body propagator %%
    % SL3 orbital propagator

    [~, oe_SL3, ~, ce_SL3] = orbprop(oe_ISSd,...
        'time',     tmax,           ...     
        'dt',       dt,             ...     
        'fmodel',   [0 0 0 0 0]     );      % No perturbation

    % Numerical integration of Kepler relative motion
    [~, oe_ODE, ce_ODE]  =  propagator01_ODE_DECHAMPS_FAYT(oe_ISSr, tspan, mu, reltol);

    % Analytical Kepler equation
    oe_ISS_KEPL = oe_ISSr;
    oe_ISS_KEPL(end) = M_ISSr;
    [~, oe_KEPL, ce_KEPL] = KEPLER_analytical(oe_ISS_KEPL, tspan, mu);
    
    % Plots comparisons
    vectors_cell = { oe_ODE; 'ODE45 integrator';
        oe_SL3; 'S3L propagator';
        oe_KEPL; 'Kepler equations'};
    keplerian_comparison(vectors_cell, tspan, MATLABc);

    % Ground tracks
    f = figure;
    f.Name = ('Ground tracks');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_ODE, 'ODE integration', 1, dt);
    subplot(2,1,2);
    grdtrk(ce_KEPL, 'KEPL propagator', 1, dt);

    orb_trk_3d(ce_ODE(:,1:3));
    final_elements_print(ce_ODE, ce_SL3, oe_ODE, oe_SL3);
    
    
elseif exo == 2

%% J2-term (Earth's oblateness) %%

    % SL3 orbital propagator
    [~, oe_SL3, ~, ce_SL3, geo_SL3] = orbprop(oe_ISSd,...
        'time',     tspan(end),  ...
        'dt',       dt,          ...
        'fmodel',   [1 0 0 0 0]);    % J2 perturbation

    
    % Numerical integration of Kepler relative motion
    [~, oe_ODE, ce_ODE, geo_ODE]  =  propagator02_ODE_DECHAMPS_FAYT(...
        oe_ISSr, tspan, mu, reltol);


    % Plots comparisons
    keplerian_comparison(oe_ODE, oe_SL3, tspan, MATLABc);
    
    % Ground track
    f = figure;
    f.Name = ('Ground tracks of the orbit');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_ODE(:,1:3), 'Fixed Earth', 0, dt);
    subplot(2,1,2);
    grdtrk(ce_ODE(:,1:3), 'Rotating Earth', 1, dt);

    
    orb_trk_3d(ce_ODE(:,1:3));
    final_elements_print(ce_ODE, ce_SL3, oe_ODE, oe_SL3);

    
elseif exo == 3
%% Earth's atmosphere %%
    % SL3 orbital propagator
    [~, oe_SL3, ~, ce_SL3] = orbprop(oe_ISSd, ...
        'time',     tspan(end),     ...
        'dt',       dt,             ...
        'fmodel',   [1 1 0 0 0],    ...     % J2 and drag perturbations
        'Cd',       ISS_Cd,         ...     
        'm',        ISS_m,          ...
        'Sd',       ISS_A,          ...
        'density',  1               );


    % Numerical integration of Kepler relative motion
    [~, oe_ODE, ce_ODE, geo_ODE]  =  propagator03_ODE_DECHAMPS_FAYT(...
        oe_ISSr, tspan, mu, ISS_prop);

    % Plots comparisons
    keplerian_comparison(oe_ODE, oe_SL3, tspan, MATLABc);

    % Ground track
    f = figure;
    f.Name = ('Ground tracks');
    f.WindowState = 'maximized';
    subplot(2,1,1);
    grdtrk(ce_ODE, 'ODE integration', 1, dt);
    subplot(2,1,2);
    grdtrk(ce_SL3, 'SL3 propagator', 1, dt);

    orb_trk_3d(ce_ODE(:,1:3));
    final_elements_print(ce_ODE, ce_SL3, oe_ODE, oe_SL3);
    
elseif exo == 4
%% Comparison with actual satellite data %%

end



%% Other functions
function final_elements_print(cart_ODE, cart_SL3, oe_ODE, oe_SL3)    
    
    % Cartesian print
    x_SL3 = cart_SL3(end,1); xdot_SL3 = cart_SL3(end,4);
    y_SL3 = cart_SL3(end,2); ydot_SL3 = cart_SL3(end,5);
    z_SL3 = cart_SL3(end,3); zdot_SL3 = cart_SL3(end,6);
    
    x_ODE = cart_ODE(end,1); xdot_ODE = cart_ODE(end,4);
    y_ODE = cart_ODE(end,2); ydot_ODE = cart_ODE(end,5);
    z_ODE = cart_ODE(end,3); zdot_ODE = cart_ODE(end,6);
    
    erel_x = abs( (x_SL3 - x_ODE) / x_SL3);
    erel_y = abs( (y_SL3 - y_ODE) / y_SL3);
    erel_z = abs( (z_SL3 - z_ODE) / z_SL3);
    
    erel_xdot = abs( (xdot_SL3 - xdot_ODE) / xdot_SL3);
    erel_ydot = abs( (ydot_SL3 - ydot_ODE) / ydot_SL3);
    erel_zdot = abs( (zdot_SL3 - zdot_ODE) / zdot_SL3);
    
    fprintf(['\nFinal cartesian coordinates are\n'             ...
       '                  SL3           Own     Error\n'       ...
       'x [km] =         %.2f     %.2f    %.1e %%\n'           ...
       'y [km] =         %.2f      %.2f     %.1e %%\n'         ...
       'z [km] =         %.2f     %.2f      %.1e %%\n'         ...
       'xdot [km/s] =    %.2f         %.2f      %.1e %%\n'     ...
       'ydot [km/s] =    %.2f         %.2f      %.1e %%\n'     ...
       'zdot [km/s] =    %.2f           %.2f      %.1e %%\n'], ...
        x_SL3/1000,    x_ODE/1000, 100 *    erel_x,            ...
        y_SL3/1000,    y_ODE/1000, 100 *    erel_y,            ...
        z_SL3/1000,    z_ODE/1000, 100 *    erel_z,            ...
     xdot_SL3/1000, xdot_ODE/1000, 100 * erel_xdot,            ...
     ydot_SL3/1000, ydot_ODE/1000, 100 * erel_ydot,            ...
     zdot_SL3/1000, zdot_ODE/1000, 100 * erel_zdot             ...
     );
    
     fprintf(['Least squared error on cartesian vectors : \n' ...
         'On position : %.2e %%\n'                            ...
         'On velocity : %.2e %%\n']                         , ...
         sqrt(erel_x^2 + erel_y^2 + erel_z^2)              , ...
         sqrt(erel_xdot^2 + erel_ydot^2 + erel_zdot^2)       ...
     );
 
    % Keplerian print
    a_SL3 = oe_SL3(end,1); w_SL3 = oe_SL3(end,4);
    e_SL3 = oe_SL3(end,2); W_SL3 = oe_SL3(end,5);
    i_SL3 = oe_SL3(end,3); t_SL3 = oe_SL3(end,6);
    
    a_ODE = oe_ODE(end,1); w_ODE = oe_ODE(end,4);
    e_ODE = oe_ODE(end,2); W_ODE = oe_ODE(end,5);
    i_ODE = oe_ODE(end,3); t_ODE = oe_ODE(end,6);
    
    erel_a = abs(a_SL3 - a_ODE) / a_SL3;
    erel_e = abs(e_SL3 - e_ODE) / e_SL3;
    erel_i = abs(i_SL3 - i_ODE) / i_SL3;
    erel_w = abs(w_SL3 - w_ODE) / w_SL3;
    erel_W = abs(W_SL3 - W_ODE) / W_SL3;
    erel_t = abs(t_SL3 - t_ODE) / t_SL3;
    
    fprintf(['\nFinal Keplerian  coordinates are\n'         ...
       '                  SL3        Own     Error\n'       ...
       'a [km] =        %.2f    %.2f    %.1e %%\n'          ...
       'e [-] =         %.2e  %.2e    %.1e %%\n'            ...
       'i [deg] =       %.2f      %.2f      %.1e %%\n'      ...
       '\x03C9 [deg] =       %.2f      %.2f      %.1e %%\n' ...
       '\x03A9 [deg] =       %.2f     %.2f     %.1e %%\n'   ...
       '\x03B8 [deg] =       %.2f     %.2f     %.1e %%\n'], ...
     a_SL3/1000, a_ODE/1000, 100 * erel_a,                  ...
          e_SL3,      e_ODE, 100 * erel_e,                  ...
          i_SL3,      i_ODE, 100 * erel_i,                  ...
          w_SL3,      w_ODE, 100 * erel_w,                  ...
          W_SL3,      W_ODE, 100 * erel_W,                  ...
          t_SL3,      t_ODE, 100 * erel_t                   ...
     );
end

function keplerian_comparison(vec_cell, tspan, MATLABc)

vec_1 = vec_cell{1};
title1 = vec_cell{2};

vec_2 = vec_cell{3};
title2 = vec_cell{4};

if numel(vec_cell) > 4
    vec_3 = vec_cell{5};
    title3 = vec_cell{6};
end

f = figure;
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
set(leg,'Position', [.455 .01 .125 .075],'Units', 'normalized');

end

function cartesian_comparison(vec_ODE, vec_SL3, tspan, MATLABc)

f = figure;
f.Name = ('Comparison of state-space vectors');
f.WindowState = 'maximized';

subplot(3,2,1)
plot( tspan/3600 ,  vec_ODE(:,1)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,1)/1000 , 'Color' , MATLABc{2}); 
title('x'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,3)
plot( tspan/3600 ,  vec_ODE(:,2)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,2)/1000 , 'Color' , MATLABc{2}); 
title('y'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,5)
plot( tspan/3600 ,  vec_ODE(:,3)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,3)/1000 , 'Color' , MATLABc{2}); 
title('z'); ylabel('Position [km]'); xlabel('Time [hours]');

subplot(3,2,2);
plot( tspan/3600 ,  vec_ODE(:,4)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,4)/1000 , 'Color' , MATLABc{2}); 
title('xdot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');

subplot(3,2,4);
plot( tspan/3600 ,  vec_ODE(:,5)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,5)/1000 , 'Color' , MATLABc{2}); 
title('ydot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');

subplot(3,2,6);
plot( tspan/3600 ,  vec_ODE(:,6)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan/3600 ,  vec_SL3(:,6)/1000 , 'Color' , MATLABc{2});
title('zdot'); ylabel('Velocity [km/s]'); xlabel('Time [hours]');


leg = legend('ode45 integrator', 'SL3 propagator');
set(leg,'Position', [.455 .01 .125 .075],'Units', 'normalized');

end


