%% Semi-major axis evolution through long propagation
% this script plots the evolution of the semi-major axis when propagation
% along 50 days is taken into account, allowing to highlight the
% atmospheric drag effect compared to the J2 perturbation.
% Each .mat file has been saved through command window after launching the
% propagators for a 50 days propagation through 1e5 time steps.

close all; clear;

MATLABc = {
    [0, 0.4470, 0.7410];
    [0.8500, 0.3250, 0.0980];
    [0.9290, 0.6940, 0.1250];
    [0.4940, 0.1840, 0.5560];
    [0.4660, 0.6740, 0.1880];
    [0.3010, 0.7450, 0.9330];
    [0.6350, 0.0780, 0.1840];
    };

load DF2_50d_1e5.mat
load DF3_50d_1e5.mat
load DF3d_50d_1e5.mat

tmax = 50*86400;
N = 1e5;
dt = tmax/N;
tspan = 0:dt:tmax;

d = 50;

f2 = figure;
plot( tspan(1:d:end)/86400 ,  oe_DF2(1:d:end,1)/1000 , 'Color' , MATLABc{1}); hold on;
plot( tspan(1:d:end)/86400 ,  oe_DF3(1:d:end,1)/1000 , 'Color' , MATLABc{2}); 
plot( tspan(1:d:end)/86400 ,  oe_DF3d(1:d:end,1)/1000, 'Color' , MATLABc{3}); 
box on; grid on; grid minor; 

ylabel('a [km]'); xlabel('Time [days]');
legend({'DF2', 'DF3', 'DF3d'}, 'Location', 'southwest');