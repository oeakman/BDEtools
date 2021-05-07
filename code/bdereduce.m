function solRed = bdereduce(sol)
% function solRed = bdereduce(sol)
%
% BDEREDUCE   Removes switch points from a BDE solution where no variables change value.
%
% solRed = bdereduce(sol)
% 
% OUTPUT
%
% solRed: A structure with spurious switch points removed.
%
% INPUTS 
%
% sol: A structure generated by e.g. bdediscret, with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%
% DEPENDENCIES
%
% None.
%
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty, University of Exeter, 2017
% k.doherty@exeter.ac.uk
%
% Code review and edits by Ozgur Akman, University of Exeter, 2019
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

%Initialise output.

solRed = struct;
solRed.x = sol.x(1);
solRed.y = sol.y(:, 1);


%Remove spurious switch points.

c = 2;

for i = 2:numel(sol.x) - 1
    
    if ~isequal(sol.y(:, i - 1), sol.y(:, i))
        
        solRed.x(c) = sol.x(1, i);
        solRed.y(:, c) = sol.y(:, i);
        
        c = c + 1;
        
    end
    
end

solRed.x = [solRed.x, sol.x(end)];
solRed.y = [solRed.y, sol.y(:, end)];

end