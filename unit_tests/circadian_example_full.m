function circadian_example_full(dataset)
% function circadian_example_full(dataset)
%
% CIRCADIAN_EXAMPLE     Demonstrates how to use BDEtools to solve a system of BDEs using a model of the clock in Neurospora crassa. The model equations are taken from 
%                       O.E. Akman et al. Digital clocks: simple Boolean models can quantitatively describe circadian systems. J. R. Soc. Interface, 9, 2365-2382 (2012).
%                       This example demonstrates a number of features of BDEtools: 
%  
%                            1. Discretising data: First, we load some synthetic data that is saved in the m-file 'neur_circ_data_gillespie.mat' (this data is the output of a
%                               run of the Gillespie algorithm) and use the bdediscrete function to convert it to logical values by thresholding.
%
%                            2. Solving with forcing input: We then use the bdesolve function to solve the equations, where we provide the light as a forcing input.
%
%                            3. Solving in serial: We also use the bdesolveserial function to solve the model where the memorisation variables are entirely given by the
%                               discretised data. Solving in serial is a useful way to check if a solution is consistent with a model. In this example, we can visually
%                               compare the output of the serial solver with the discretised data and observe that they are similar. Hence, we can conclude that the data is
%                               consistent with the model. This could be a particularly useful test if our parallel prediction did not fit well with our data.
% 
%                            4. Manipulating and plotting solutions: We use the bdejoin function to combine the history with a prediction into one solution. We use the
%                               bdemerge function to combine the light forcing with the model prediction into one solution. We use the bdePR function to convert the solutions to
%                               straight line segments, to facilitate easy plotting.
%
% circadian_example_full('deterministic')
% circadian_example_full('stochastic')
%
% OUTPUT
%
% Plots of the discretised data and the parallel/serial solutions.
%
% INPUTS
%
% dataset: Chooses whether a deterministic or stochastic dataset is used.
%
% DEPENDENCIES 
%
% bdediscrete, bdesolve, bdesolveserial, bdejoin, bdemerge, bdePR, bdedist.
%
% SEE ALSO
%
% bdesolve, bdesolveserial.
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty & Ozgur Akman, University of Exeter, 2021
% k.doherty@exeter.ac.uk
% O.E.Akman@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2021
%

% Specify final timepoint of prediction.

tEnd = 120;

% Define the model parameters (delays). The following are the average times between switches in the data (these have been calculated in advance).

tau1 = 5.0752; 
tau2 = 6.0211;
tau3 = 14.5586;
lags = [tau1, tau2, tau3];

% Load data to be used for the history.

if strcmp(dataset, 'stochastic')
    
    load('neur_circ_data_gillespie.mat'); % Load the synthetic dataset generated using the Gillespie algorithm.
    
    xData = neur_circ_data_gillespie.LD(3, :);
    yData = neur_circ_data_gillespie.LD(1:2, :);
    
elseif strcmp(dataset, 'deterministic')
    
    load('neur_circ_data_ode.mat'); % Load the synthetic dataset generated from the equivalent set of ODEs.
    
    xData = neur_circ_data_ode.LD(3, :); 
    yData = neur_circ_data_ode.LD(1:2, :);
    
end

T = [0.3 0.3]; % Set thresholds for discretising the data.

synData = bdediscrete(xData, yData, T); % Convert real-valued data to Boolean data by thresholding.

if synData.x<tEnd
    synData.x = [synData.x, tEnd]; % Add an extra point at the end to give the data the same tRange as the prediction.
    synData.y = [synData.y, synData.y(:, end)]; 
end

% Specify the history.

history.x = synData.x(synData.x < 24); % Define the history as all points from t=0 to t=24.
history.y = synData.y(:, synData.x < 24);

% Augment the history.

history.x = [history.x, 24]; % Include a point at t=24 to make the history a full 24 hours (note: only the values from 24 - max(lags) to 24 will be used in calculating the solution).
history.y = [history.y, history.y(:, end)];

% Specify the forcing.

forcing.x = [0, 6, 18:12:114, 120]; 
forcing.y = mod(forcing.x, 24) >= 6 & (mod(forcing.x, 24) < 18);

% Plot the data.

figure;
plot_neurospora_data(xData, yData, T, synData.x, synData.y);

% Solve the equations in parallel.

solPar = bdesolve(@neurospora_eqns, lags, history, tEnd, forcing); 

% Solve the equations in serial.

solDataSer = bdesolveserial(@neurospora_eqns, lags, synData, [24 tEnd], forcing);

% Append the histories.

solDataSer_wHist = bdejoin(solDataSer.history, solDataSer);
solPar_wHist = bdejoin(solPar.history, solPar);

% Solve the equations in serial, with the parallel solution as the history.

solSer = bdesolveserial(@neurospora_eqns, lags, solPar_wHist, [24 tEnd], forcing);
solSer_wHist = bdejoin(solSer.history, solSer);

% Calculate Hamming distances.

d1 = sum(bdedist(solDataSer_wHist, synData)/tEnd); 
d2 = sum(bdedist(solPar_wHist, synData)/tEnd);
d3 = sum(bdedist(solSer_wHist, solPar_wHist)/tEnd);

disp(sprintf('Normalised Hamming distance between serial solution and the data: \n %d', d1)); %#ok
disp(sprintf('Normalised Hamming distance between parallel solution and the data: \n %d', d2)); %#ok
disp(sprintf('Verify solution: Normalised Hamming distance between the parallel \n solution and the solution obtained by running it through the serial solver: \n %d', d3)); %#ok

% Plot the simulations.

figure;
subplot(2, 1, 1);
plot_neurospora_prediction(solPar);
title('Model prediction (parallel)');
subplot(2, 1, 2);
plot_neurospora_prediction(solDataSer)
title('Model prediction (serial)');

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction that defines the model.

function X = neurospora_eqns(Z1, Z2)

M = ~Z1(2, 2) | Z2(1, 3);
FT = Z1(1, 1);

X = [M; FT];

% Subfunction for plotting the continuous and discretised data.

function plot_neurospora_data(xData, yData, T, xDataBool, yDataBool)

offset = 1.1;

data.x = xDataBool;
data.y = yDataBool;

dataPR = bdePR(data, offset);

T1 = min(yData(1, :)) + T(1) * (max(yData(1, :)) - min(yData(1, :)));
T2 = min(yData(2, :)) + T(2) * (max(yData(2, :)) - min(yData(2, :)));

my_colours = get(gca,'colororder');

subplot(2, 1, 1);
plot(xData, yData, 'Marker', 'o', 'LineStyle', 'none', 'LineWidth', 2); % Plot the data. 
legend('FRQ mRNA', 'FRQ Protein','FRQ threshold','FRQ Protein');
hold on ;
plot([0 120], [T1, T1], 'Color', my_colours(1, :), 'LineWidth', 2);
plot([0 120], [T2, T2], 'Color', my_colours(2, :), 'LineWidth', 2);
xlabel('Time, t (h)');
ylabel('Concentration (nM)');
title('Neurospora crassa circadian data (with thresholds)');

subplot(2, 1, 2);
plot(dataPR.x, dataPR.y, 'LineWidth', 2) % Plot discretised data.
hold on;
plot([24 24], [0 3.3], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2) % Also plot the line separating the history from the prediction.
for i = 24:24:120 % And plot lines to separate days.
    plot([i i], [0 3.3], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
end
yticks([0.5, 1.6]); % Have one ytick for each variable.
yticklabels({'F_M(t)', 'F_P(t)'}); % Name each ytick.
ylim([0 2.2]);
xlabel('Time, t (h)');
ylabel('Activity (1/0)');
title('Thresholded data');

% Subfunction for plotting model predictions. 

function plot_neurospora_prediction(sol)

solPlot = bdejoin(sol.history, sol);
solPlot = bdemerge(solPlot, sol.forcing);

offset = 1.1;
tEndHist = sol.history.x(end);
solPR = bdePR(solPlot, offset);

plot(solPR.x, solPR.y, 'LineWidth', 2);
hold on;

for i = 24:24:120 % And plot lines to separate days.
    plot([i i], [0 3.3], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
end
plot([tEndHist tEndHist], [0 3.3], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2); % Also plot the line separating the history from the prediction.
yticks([0.5, 1.6, 2.7]);
yticklabels({'F_M(t)','F_T(t)','L(t)'});
ylim([0 3.3]);
xlabel('Time t, (h)');
ylabel('Activity (1/0)');






