function sol = ENSO_example(tEnd)
%function sol = ENSO_example(tEnd)
%
% ENSO_EXAMPLE     A Boolean Delay Equation model of the El Nino/Southern Oscillation phenomenona. The equations are taken from 
%                  A. Saunders & M. Ghil. A Boolean delay equation model of ENSO variability. Physica D, 160:54?78 (2001),
%                  and the systen is solved for the parameter values tau = 0.31, beta = 0.77 and seasonal forcing period, P = 2. 
%                  This replicates Fig. 2(a)i). Model variables:
%                       U1 represents westerly wind anomolies.
%                       U2 represents easterly wind anomolies.
%                       T1 represents cold/hot (0/1) SSTA anomolies.
%                       T2 represents amplification (if equal to T1).
%                       S represents seasonal forcing. 
%                  The dynamics are summarised using the variable, ENSO, with four distinct states:
%                       E = -2: Extremely cold SSTA (La Nina)
%                       E = -1: Mild cold SSTA
%                       E = 1: Mild warm SSTA
%                       E = 2: Extremely warm SSTA (El Nino)
% 
% sol = ENSO_example(tEnd)
%
% OUTPUT
%
% sol: A structure with the following fields - 
% sol.x: A vector containing the times of switch points.
% sol.y: A matrix with 2 rows where the rows are the state variables. Each column is the state following each switch.
% sol.solver: A string specifying the solver used (bdesolve).
% sol.bdefun: The system of equations that is solved.
% sol.lags: A vector of lag values.
% sol.numEvals: The number of times the equations were called.
% sol.history: The history used to generate the solution.
% sol.historyFlag: A logical value which specifies whether a switch occurs at the first point in the solution. This can indicate if the 
%                  history is incompatible with the solution. 
% sol.earlyTermFlag: A logical value which specifies whether the solution was terminated prematurely (e.g. if the maximum number of 
%                    switch points was exceeded). 
%
% INPUTS
%
% tEnd: Sets the timespan - the equations are integrated over [0 tEnd].
%
% DEPENDENCIES 
%
% bdesolve, bdePR.
%
% SEE ALSO
%
% bdesolve.
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty & Ozgur Akman, University of Exeter, 2021
% k.doherty@exeter.ac.uk
% O.E.Akman@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2021
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

% Set the parameter values.

tau = 0.31;
beta = 0.17;
P = 2;

% Set the delays.

lags = [tau, beta, P]; 

% The following history gives the limit cycle solution shown in the paper 
% (this is sensitive to the initial conditions).

history = struct;
history.x = [0 0.15 0.32 0.83 1 1.17 1.48 1.5 1.65 1.67 1.84 2];
U1h = [false, false, true, true, true, true, true, true, false, false, false, false];
U2h = [false, false, false, false, false, true, true, true, true, true, false, false];
T1h = [false, true, true, true, true, true, false, false, false, false, false, false];
T2h = [false, false, false, false, true, true, true, true, true, false, false, false];
Sh = [false, false, false, true, true, true, true, false, false, false, false, false];
history.y = [U1h; U2h; T1h; T2h; Sh];

% Solve the equations.

sol = bdesolve(@ENSO_eqns, lags, history, tEnd);

% Plot the solution.

plot_ENSO(sol);

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction to evaluate the BDEs at time t given the memorisation variable values. 
% 
% The matrix Z contains the memorisation variables, i.e. the values of the model variables at the current time t minus each of the lags.
% 
% Z = [U1(t - tau), U1(t - beta), U1(t - P);
% U2(t - tau), U2(t - beta), U2(t - P);
% T1(t - tau), T1(t - beta), T1(t - P); 
% T2(t - tau), T2(t - beta), T2(t - P)];
%

function X = ENSO_eqns(Z)

U1 = Z(1,:); 
U2 = Z(2,:); 
T1 = Z(3,:); 
T2 = Z(4,:);
S = Z(5,:);

Delta = @(x1, x2) ~xor(x1, x2); % Returns TRUE if both inputs equal.

R = Delta(U1(1), U2(1)); % Rossby-wave signal.

U1t = T1(2); % U1(t) = T1(t - beta).

U2t = T2(2); % U2(t) = T2(t - beta).

T1t = (R & ~U1(1)) | (~R &  U1(2)); 

T2t = S(2); % Periodicity is determined by seasonal forcing.

St = S(3); % S(t) = S(t - P). Seasonal forcing term. 

X = [U1t; U2t; T1t; T2t; St];

end

% Subfunction to plot the history and the model prediction. A vertical dashed line separates the history from the prediction. 

function plot_ENSO(sol)

% Get the point where history meets the prediction.
    
histEnd = sol.history.x(end); 

% Offset solutions for plotting.

offset = 1.1; 

% Join the history and prediction into one timeseries solution.

tSeries1 = bdejoin(sol.history, sol); 

% Convert to lines for plotting.

tSeries = bdePR(tSeries1); 

% Calculate the ENSO variable.

ENSO = zeros(1, numel(tSeries.x));
ENSO(tSeries.y(3, :) == 1 & tSeries.y(4, :) == 0) = 1; % If T1 == 1 then ENSO positive.
ENSO(tSeries.y(3, :) == 1 & tSeries.y(4, :) == 1) = 2;
ENSO(tSeries.y(3, :) == 0 & tSeries.y(4, :) == 0) = -2; % If T1 == 0 then ENSO negative.
ENSO(tSeries.y(3, :) == 0 & tSeries.y(4, :) == 1) = -1;

% Plot the ENSO values.

tSeries = bdePR(tSeries1, offset);
figure;
subplot(2, 1, 1);
hold on;
plot(tSeries.x, ENSO,'LineWidth',2);
plot([histEnd histEnd], [-2.2 2.2], '--', 'LineWidth', 2, 'Color', 'k');
xlabel('Time, t (years)');
xlim([tSeries.x(1) tSeries.x(end)]);
ylim([-2.2 2.2]);
yticks([-2 -1 0 1 2]);
title('ENSO');
box on;

% Plot each variable separately also.

subplot(2, 1, 2);
hold on;
plot(tSeries.x, tSeries.y, 'LineWidth', 2);
plot([histEnd histEnd], [0 5*offset], '--', 'LineWidth', 2, 'Color', 'k');
title('Individual variables');
xlim([tSeries.x(1) tSeries.x(end)])
ylim([0 5*offset]);
xlabel('Time, t (years)');
yticks(offset.*[0.5, 1.5, 2.5, 3.5, 4.5]);
yticklabels({'U1','U2','T1','T2','S'});
box on;

end

end