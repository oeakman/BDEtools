function [sol, solSer, solerr] = chaotic_mod_test_wplot
% function [sol, solSer, solerr] = chaotic_mod_test_wplot
%
% CHAOTIC_MOD_TEST_WPLOT   Test function for the BDE solver - chaotic example of Dee et al. 
%                          This version also plots the parallel & serial solutions for comparison.
%
% [sol, solSer, solerr] = chaotic_mod_test
% [~, ~, solerr] = chaotic_mod_test
%
% OUTPUTS 
%
% sol/solSer: Structures returned by bdesolve/bdesolveserial, with the following common fields -
% sol(Ser).x: A vector containing the times of switch points.
% sol(Ser).y: A matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% sol(Ser).history: The history used to generate the solution.
% 
% solerr: Error between the parallel and serial solutions (calculated using continuous normalised Hamming distance).
%
% INPUTS
%
% None.
%
% MODEL
%
% The model equations are:
%
% x_1(t) = x_2(t - tau_1),
% x_2(t) = x_1(t - tau_1) XOR x_2(t - tau_2).
%
% Reference: Dee, D. & Ghil, M. SIAM J. Appl. Math. 44(1): 111?126 (1984).
%
% DEPENDENCIES 
%
% bdesolve.
% bdesolveserial.
% bdecut.
% bdedist.
% bdeplot.
% 
% SEE ALSO 
%
% bdesolve.
% bdesolveserial.
% chaotic_mod_test.
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty & Ozgur Akman, University of Exeter, 2017
% k.doherty@exeter.ac.uk
% O.E.Akman@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2019
%

% Set the delays.

tau1 = 0.977;
tau2 = 1;

% Set the range over which to integrate.

tRange = [0 50];

% Set the (constant) history.

history = [true; true];

% Define the delays.

lags = [tau1 tau2];

% Specify the model to use.

fun = @chaotic_eqns;

% Solve the equations.

sol = bdesolve(fun, lags, history, tRange);

% Solve the equations serially.

solSer = bdesolveserial(fun, lags, sol, tRange + max(lags));

% Extract the coincident regions of the two solutions.

[~, solCut] = bdecut(sol, solSer.x(1));
[solSerCut, ~] = bdecut(solSer, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr] = bdedist(solCut, solSerCut);

% Plot the parallel and serial solutions.

figure; 
subplot(2, 1, 1);
bdeplot(sol, 1.1);
subplot(2, 1, 2);
bdeplot(solSerCut, 1.1);

for k=1:2
    subplot(2, 1, k);
    ax = axis;
    ax(3:4) = [-0.05 2.5];
    axis(ax);
    box on;
    xticks(0:5:50); 
    yticks([0 1 1.1 2.1]); 
    yticklabels({'0', '1', '0', '1'}); 
    set(gca, 'FontSize', 12);
    legend({'x_1', 'x_2'}, 'Orientation', 'Horizontal', 'Location', 'North', 'FontSize', 18);
end

subplot(2, 1, 1);
title('Parallel solution', 'FontSize', 20);
subplot(2, 1, 2);
title('Serial solution', 'FontSize', 20);

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction defining the model.

function X = chaotic_eqns(Z)

X = [Z(2, 1);
    xor(Z(1, 1), Z(2, 2))];

end

end

