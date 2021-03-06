function [per, solper] = bdeper(sol)
% function [per, solper] = bdeper(sol)
%
% BDEPER   Computes a single period of a limit cycle BDE solution.
%
% per = bdeper(sol)
% [per, solper] = bdeper(sol)
% 
% OUTPUTS
%
% per: Approximate period of the solution.
% solper: Structure containing one period of the solution (OPTIONAL).
% 
% INPUTS 
%
% sol: A structure generated by e.g. bdesolve, with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%
% DEPENDENCIES
%
% None.
%
% -------------------------------------------------------------------------
%
% Written Ozgur Akman, University of Exeter, 2019
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

% Convert the solution into a character string.

ymat = double(sol.y');
ymat = num2str(ymat);

% Convert the bit strings into decimals and look for periodicity.

bstr = bin2dec(ymat);
in1 = (sum(bstr == bstr'));
in2 = find(in1 == min(in1));
in = find(bstr == bstr(in2(1)));

if length(in) > 2
   
    per = sol.x(in(end)) - sol.x(in(end-1));
    
    % If requested, return a single period of the oscillation.
    
    if nargout > 1
        
        solper = sol;
        solper.x = sol.x(in(end-1):in(end));
        solper.x = solper.x - solper.x(1);
        solper.y = sol.y(:,in(end-1):in(end));       
        
    end
        
else
    
    per = Inf;
       
end

end