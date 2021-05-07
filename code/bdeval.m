function y = bdeval(sol, x)
% function y = bdeval(sol, x)
%
% BDEVAL    Returns the values of a BDE solution at some specfied timepoints. 
% 
% y = bdeval(sol, x)
%
% OUTPUT
%
% y: The values of the BDE solutions at the required points.
% 
% INPUTS
%
% sol: A structure with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n1 rows where n1 is the number of state variables. Each column is the state following each switch.
% 
% DEPENDENCIES 
%
% None.
%
% SEE ALSO
%
% bdesolve.
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty & Ozgur Akman, University of Exeter, 2017
% k.doherty@exeter.ac.uk
% O.E.Akman@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2020
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

if any(x < sol.x(1)) || any(x > sol.x(end))    
    error('The values of x must be within the range of sol.x.')    
end

idxLT = sol.x <= x(:);
y = sol.y(:, sum(idxLT, 2));

end

