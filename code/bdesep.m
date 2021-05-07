function solSep = bdesep(sol, idx)
% function solSep = bdesep(sol, idx)
%
% BDESEP    Extract variables from a BDE solution.
% 
% solSep = bdesep(sol, idx)
%
% OUTPUT
%
% solSep: A structure with the following fields - 
% solSep.x: A vector with the times of switch points for the variables of interest.
% solSep.y: A Boolean matrix with m rows, where m is the number of variables of interest. Each column is the state following each switch.
% 
% INPUTS
%
% sol: A structure with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n rows where n (>=m) is the number of state variables. Each column is the state following each switch.
% 
% idx: A vector of integers indexing the variables to be extracted.
%
% DEPENDENCIES 
%
% None.
%
% SEE ALSO
%
% bdemerge (inverse function).
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty & Ozgur Akman, University of Exeter, 2017
% k.doherty@exeter.ac.uk
% O.E.Akman@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2019
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

% Get indices of any columns of sol.y where there is no change compared with the previous column (doesn't include last point).

bla = any(xor(sol.y(idx, 2:end-1), sol.y(idx, 1:end-2)),1); 

% Extract variables.

solSep.x = sol.x([true, bla, true]);
solSep.y = sol.y(idx, [true, bla, true]);

end
