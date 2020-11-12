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

%% Convergence study

% Relative tolerances array
reltol = ones(1,10) * 1e-3;
for i = 1:numel(reltol)
    reltol(i) = reltol(i)*10^(-i);
end

% Norms of cartesian relative errors (computed separately)
norm_erel_cart = [
       2.49     1.41 ;
    8.38E-2  3.62E-2 ;
    2.36E-4  1.01E-4 ;
    1.82E-4  7.86E-5 ;
    2.35E-5  1.02E-5 ;
    3.11E-6  1.35E-6 ;
    9.38E-7  4.25E-7 ;
    7.72E-7  3.34E-7 ;
    7.51E-7  3.25E-7 ;
    7.49E-7  3.24E-7 ]; 

% Plot
figure; box on; grid on;

plot(reltol,norm_erel_cart(:,1),  ...
    'Color', MATLABc{1},          ...
    'LineWidth', 1 );
hold on;
plot(reltol,norm_erel_cart(:,2),  ...
    'Color', MATLABc{2},          ...
    'LineWidth', 1,               ...
    'LineStyle', '--'              );

set(gca,                ...
    'XDir', 'reverse',  ...
    'YScale', 'log',    ...
    'XScale', 'log'      );

xlabel('Relative tolerance'); 
ylabel('Norm of relative error vector');
xlim([1e-15 1e-5]);
legend({'Position', 'Velocity'}),

