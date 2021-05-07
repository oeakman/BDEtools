function sol = bdemerge(sol1, sol2)
% function sol = bdemerge(sol1, sol2)
%
% BDEMERGE    Merge two BDE solutions defined over the same range.
% 
% sol = bdemerge(sol1, sol2)
%
% OUTPUT
%
% sol: A structure with the following fields - 
% sol.x: A vector with the times of combined switch points for the variables in the two solutions.
% sol.y: A Boolean matrix with m (=n1+n2) rows, where m is the number of combined variables in the two solutions. Each column is the state following each switch.
% 
% INPUTS
%
% sol1: A structure with the following fields - 
% sol1.x: A vector with the times of switch points.
% sol1.y: A Boolean matrix with n1 rows where n1 is the number of state variables. Each column is the state following each switch.
% 
% sol2: A structure with the following fields - 
% sol2.x: A vector with the times of switch points.
% sol2.y: A Boolean matrix with n2 rows where n2 is the number of state variables. Each column is the state following each switch.
%
% DEPENDENCIES 
%
% None.
%
% SEE ALSO
%
% bdesep (inverse function).
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

if numel(fieldnames(sol1)) == 0
    
    sol = sol2;
    
    return
    
elseif numel(fieldnames(sol2)) == 0
    
    sol = sol1;
    
    return
    
end

if (sol1.x(1) ~= sol2.x(1)) || (sol1.x(end) ~= sol2.x(end))
    
    error('The two solutions (sol1, sol2) must have the same range for their independent variables (x).')
    
end

sol.x = []; 
sol.y = logical([]); 

prevIdxGt = 1;

for i = 2:numel(sol2.x)
    
    idxGt = find(sol1.x > sol2.x(i));
        
    if isempty(idxGt)
        
        r = prevIdxGt : numel(sol1.x);
        
    else
        
        idxGt1 = idxGt(1);
        r = prevIdxGt : idxGt1 - 1;
        
        prevIdxGt = idxGt1;
        
    end
    
    if isempty(r)
        
        sol.x = [sol.x, sol2.x(i)];
        sol.y = [sol.y, [sol1.y(:, idxGt1 - 1); sol2.y(:, i)]];
        
    elseif isequal(sol1.x(r(end)), sol2.x(i))
        
        sol.x = [sol.x, sol1.x(r)];
        sol.y = [sol.y, [sol1.y(:, r); [repmat(sol2.y(:, i-1), 1, numel(r)-1), sol2.y(:, i)]]];
        
    else
        
        sol.x = [sol.x, sol1.x(r), sol2.x(i)];
        sol.y = [sol.y, [sol1.y(:, r); repmat(sol2.y(:, i-1), 1, numel(r))], [sol1.y(:, idxGt1 - 1); sol2.y(:, i)]];
        
    end
    
end

end

