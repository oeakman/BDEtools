function sol = chaotic_example(tEnd)
%function sol = chaotic_example(tEnd)
%
% CHAOTIC_EXAMPLE     The simple Boolean Delay Equation model of Dee & Ghil. The equations are taken from 
%                     D. Dee & M. Ghil. Boolean difference equations, i: Formulation and dynamic behavior. SIAM J. Appl. Math., 44(1):111?126 (1984),
%                     and the systen is solved for the parameter values tau_1 = 0.977, tau_2 = 1.
% 
% sol = chaotic_example(tEnd)
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

% Set the delays.

tau1 = 0.977;
tau2 = 1;

% Set the timespan.

tRange = [0 tEnd];

% Set the history.

history = [true; true];

% Define the BDE functions.

lags = [tau1 tau2];
fun = @chaotic_eqns;

% Solve the equations. 

sol = bdesolve(fun, lags, history, tRange);

% Plot the solution.

solPR = bdePR(sol, 1.1);

figure;
plot(solPR.x, solPR.y, 'LineWidth', 2);
yticks([0.5, 1.6]); 
yticklabels({'x_1', 'x_2'});
xlabel('Time, t');

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction defining the equations. 

function X = chaotic_eqns(Z)

X = [Z(2, 1);
    xor(Z(1, 1), Z(2, 2))];

end

end

