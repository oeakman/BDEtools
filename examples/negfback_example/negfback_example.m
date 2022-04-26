% NEGFBACK_EXAMPLE     Compares the analytical and numerical solutions to a Boolean Delay Equation model of a simple negative feedback circuit.
% 
% negfback_example
%
% OUTPUT
%
% None.
%
% INPUTS
%
% None.
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

% Generate and plot the analytical solution.

% Specify the time points of switches

asol.x = [0 1 2 3 4 5 5.5]; 			

% Specify the states after each switch in asol.x

asol.y = [[0;0], [1;0], [1;1], [0;1], [0;0], [1;0], [1;0]]; 	

% Make the solution plot-ready 
% (i.e. convert to straight line segments)

asolPR = bdePR(asol, 1.1); 						
figure;
plot(asolPR.x, asolPR.y, 'LineWidth', 2);
xlim([0 5.5]);
ylim([0 2.3]);
yticks([0.5 1.6]);
yticklabels({'A', 'B'});
xlabel('Time, t');

% Generate and plot the numerical solution for comparison.

tau1 = 1;
tau2 = 1;

% Define the delays

lags = [tau1, tau2]; 

% Constant history for all time t<0

history = [false; true]; 

% Generate solution for 5.5 units of time

tspan = [0 5.5]; 

% The model equations, where Z contains the values of the 
% memorisation variables.

fun = @(Z) [~Z(2,2); Z(1,1)]; 
nsol = bdesolve(fun, lags, history, tspan);

% Make the solution plot-ready

nsolPR = bdePR(nsol, 1.1); 
figure;
plot(nsolPR.x, nsolPR.y, 'LineWidth', 2);
xlim([0 5.5]);
ylim([0 2.3]);
yticks([0.5 1.6]);
yticklabels({'A', 'B'});
xlabel('Time, t');