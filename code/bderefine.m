function solRef = bderefine(sol, yvar, xtol)
% function solRef = bderefine(sol, yvar, xtol)
%
% BDEREFINE   Removes successive switch points from a selected variable of a BDE solution, when the switches occur faster than a specified tolerance. 
%
% solRef = bderefine(sol, yvar, xtol)
% solRef = bderefine(sol, yvar, []) 
% solRef = bderefine(sol, yvar)
%
% OUTPUT
%
% solRef: A structure with the switch points removed.
%
% INPUTS 
%
% sol: A structure generated by, e.g., bdediscrete, with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%
% yvar: An integer indexing the selected variable. 
% xtol: The specified tolerance (OPTIONAL: default value is 1e-6).
%
% DEPENDENCIES
%
% bdesep, bdemerge.
%
% -------------------------------------------------------------------------
%
% Written by Ozgur Akman, University of Exeter, 2022
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

% Set the default value of the tolerance.

if (nargin < 3) || (isempty(xtol))
    xtol = 1e-6;
end

%Extract the required variable.

sol1 = bdesep(sol, yvar);

%Extract the remaining variables.

sol2 = bdesep(sol, setdiff(1:size(sol.y,1), yvar));

%Remove switch points in which the time between successive switches is less than the tolerance.

dx = diff(sol1.x);
in = find(dx<xtol);

%If any of these occur at the beginning or end of the time interval, ignore.

in(in<2) = [];
in(in>(length(sol1.x)-2)) = [];

if ~isempty(in)
    din = [in in+1];
    sol1.x(din) = [];
    sol1.y(:, din) = [];
end

%Merge the solutions.

solRef = bdemerge(sol1, sol2);
solRef.xtol = xtol;

end