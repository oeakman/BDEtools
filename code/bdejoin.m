function solJ = bdejoin(sol1, sol2)
% function solJ = bdejoin(sol1, sol2)
%
% BDEJOIN   Splices together two BDE solutions in which the end point of one solution is the initial point of the other, 
%           and which share the same # of variables.
%
% solJ = bdejoin(sol1, sol2)
%
% OUTPUT
%
% solJ: A structure containing the spliced solution, with the following fields - 
% solJ.x: A vector with the times of switch points.
% solJ.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
 %
% INPUTS 
%
% sol1, sol2: Structures with the following fields - 
% sol[1/2].x: A vector with the times of switch points, with sol1.x(end)=sol2.x(1).
% sol[1/2].y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%  
% DEPENDENCIES
%
% None.
%
% SEE ALSO.
%
% bdecut (inverse function).
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

% Check that sol1 and sol2 are compatible.

if sol1.x(end) ~= sol2.x(1)
    
    error('The final value of sol1.x should equal the first value of sol2.x.')
    
end

if size(sol1.y, 1) ~= size(sol2.y, 1)
    
    error('sol1 and sol2 should have the same number of rows.')
    
end

% Splice the solutions together.
   
if sum(xor(sol2.y(:, 1), sol2.y(:, 2)))% If there is a switch at the first point in sol2, we need to include this point.
    
    solJ.x = [sol1.x(1: end - 1), sol2.x];
    solJ.y = [sol1.y(:, 1: end - 1), sol2.y];
    % [Note that a switch in the second solution takes priority, i.e. will overwrite a switch if it occurs at the end of sol1.]
    
elseif sum(xor(sol1.y(:, end), sol1.y(:, end - 1))) % If there is a switch at the last point in sol1 we need to include this point.
    
    solJ.x = [sol1.x, sol2.x(2:end)];
    solJ.y = [sol1.y, sol2.y(:, 2:end)];
    
else
    % Otherwise, we remove the point that joins them.
    
    solJ.x = [sol1.x(1:end - 1), sol2.x(2:end)];
    solJ.y = [sol1.y(:, 1:end - 1), sol2.y(:, 2:end)];
    
end

end