function [sol, solSer, solerr] = chaotic_mod_test
% function [sol, solSer, solerr] = chaotic_mod_test
%
% CHAOTIC_MOD_TEST   Test function for the BDE solver - chaotic example of Dee et al.
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

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction defining the model.

function X = chaotic_eqns(Z)

X = [Z(2, 1);
    xor(Z(1, 1), Z(2, 2))];

end

end

