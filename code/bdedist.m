function [d, nd, mnd] = bdedist(sol1, varargin)
% function [d, nd, mnd] = bdedist(sol1, varargin)
% 
% BDEDIST   Calculate continuous Hamming distance.
% 
% d = BDEDIST(sol1)
% d = BDEDIST(sol1, sol2) 
% [d, nd] = BDEDIST(sol1, sol2) 
% [d, nd, mnd] = BDEDIST(sol1, sol2) 
%
% OUTPUTS
%
% d: Vector containing either the Hamming distances of each variable in sol1 from 0 or the distances between the variables in sol1 and 
%    sol2 (entered as varargin{1} - see below).
% 
% nd: The normalised Hamming distances for each variable (OPTIONAL).
%
% mnd: The mean normalised Hamming distance across all variables (OPTIONAL).
%
% INPUTS
%
% sol1: A structure with the following fields - 
% sol1.x: A vector with the times of switch points.
% sol1.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% 
% varargin{1} -
% sol2: A structure with the following fields - 
% sol2.x: A vector with the times of switch points.
% sol2.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%
% DEPENDENCIES
%
% None.
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

% Case where required distance is between two solutions.

if nargin > 1
    
    sol2 = varargin{1};
    
    if (sol1.x(1) ~= sol2.x(1)) || (sol1.x(end) ~= sol2.x(end))
        
        error('Both solutions should be defined over the same range of values of x.')
        
    end
    
else
    
    sol2.x = [sol1.x(1) sol1.x(end)];
    sol2.y = false(size(sol1.y, 1), 2);
    
end

[switchPointsSorted, idxSort] = sort([sol1.x sol2.x]); % Sort all switches (from both solutions).

% Define a set of indices to define if a switch is from sol1 or sol2.

solIdx = [false(1, numel(sol1.x)), true(1, numel(sol2.x))]; 
solIdx = solIdx(idxSort); 

segs = diff(switchPointsSorted); % Get length of each segment.

d = zeros(size(sol1.y, 1), 1);

c1 = 1; % Counters for the switches in sol1 and sol2, respectively.
c2 = 1;

for i = 1:numel(segs)-2 % Step through segments.
    
    d = d + xor(sol1.y(:, c1), sol2.y(:, c2)) * segs(i+1); % Add the length of each segment to the total distance where solutions differ.
    
    if ~solIdx(i+2) % If a switch occurs in sol1.
        
        c1 = c1 + 1;
        
    else % If a switch occurs in sol2.
        
        c2 = c2 + 1;
        
    end
            
end

% Calculate the normalised & mean Normalised Hamming distances, if required.

if nargout>1
    
    nd = d / (sol1.x(end) - sol1.x(1));
    
    if nargout > 2
       
        mnd = mean (nd);
        
    end

end
  
end
       