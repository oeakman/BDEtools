function sol = bdesolveserial(bdefun, lagsRVec, serialData, xRange, varargin)
% function sol = bdesolveserial(bdefun, lagsRVec, serialData, xRange, varargin)
%
% BDESOLVESERIAL   Solve a system of Boolean Delay Equations in series. Works similarly to bdesolve.m 
%                  but memorisation variable values are taken directly from serialData. In other words, 
%                  the state of each variable in the solution depends only on the input data.
% 
% sol = BDESOLVESERIAL(bdefun, lagsRVec, serialData, xRange)
% sol = BDESOLVESERIAL(bdefun, lagsRVec, serialData, xRange, forcing)
%
% OUTPUT 
%
% sol: A structure with the following fields -
% sol.x: A vector containing the times of switch points.
% sol.y: A matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% sol.history: The portion of <SerialData> preceding the projection interval.
% sol.lagsRvec: A vector of lag values.
% sol.solver: A string specifying the solver used.
% sol.bdefun: The system of equations that is solved.
% sol.forcing: This is only present if forcing was provided.
%
% INPUTS
%
% bdefun: Handle to function that takes in a matrix of memorisation variables (each column corresponds to each value in <lags>) and returns 
%         a vector. If there is forcing, the function also takes in a matrix of memorisation variables for the forcing (see below).
%
% lagsRvec: A vector of lag values.
%
% serialData: A structure that should at least contain the following two fields -
% serialData.x: A vector containing the times of switch points.
% serialData.y: A matrix with n rows where n is the number of state variables. Each column is the state following each switch.
%
% xRange: A 1 x 2 vector containing the start and end points of integration.
%
% varargin{1} -
% forcing: A structure with two fields, forcing.x and forcing.y. forcing.x should span the same range of values as xRange. If forcing is 
% provided, then bdefun should take an extra argument, the memorisation variables for the forcing.
%
% DEPENDENCIES
% 
% bdecut.
%
% SEE ALSO 
%
% bdesolve (the parallel solver).
% 
% -------------------------------------------------------------------------
%
% Written by Kevin Doherty, University of Exeter, 2017
% k.doherty@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2019
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

potSwitches = [];
memDelaysIdx = [];
memVarsIdx = [];
numLags = numel(lagsRVec);

% Initialise forcing.

if nargin > 4
    
    if ~isempty(varargin{1})
        
        INCL_FORCING = true;
        forcing = varargin{1};
        
        numForceSwitches = numel(forcing.x);
        sol.forcing = forcing;

        for idxForceSwitches = 1:numForceSwitches
        
            potSwitches = [potSwitches, lagsRVec + forcing.x(idxForceSwitches)];
        
        end
        
        memVarsIdx = zeros(1, numForceSwitches * numLags); % Includes a value of 0 for any switches caused by forcing.
        memDelaysIdx = zeros(1, numForceSwitches * numLags); 
    
        
    else
        
        INCL_FORCING = false;
        
    end
    
else
    
    INCL_FORCING = false;
    
end

% ----------- %

for idxSwitch = 1:numel(serialData.x)
    
    potSwitches = [potSwitches, serialData.x(idxSwitch) + lagsRVec];
    memDelaysIdx = [memDelaysIdx, 1:numLags];
    memVarsIdx = [memVarsIdx, idxSwitch * ones(1, numLags)];
    
end

xEnd = xRange(end);

if numel(xRange) > 1
    
    x0 = xRange(1);
    
    if x0 < serialData.x(1) + max(lagsRVec)
        
        error('The initial point in the range does not have a defined history.')
        
    end
   
else
    
    x0 = max(lagsRVec);
    
end


[sol1, sol2] = bdecut(serialData, x0);

sol.x = sol1.x;
sol.y = sol1.y;

lenHist = numel(sol.x) - 1;

sol.history = sol;
sol.solver = 'bdesolveserial';
sol.bdefun = bdefun;

xNow = sol2.x(1);

someIdx = find(potSwitches == xNow);
    
% --------------------------- %
% ----- Begin iteration ----- %
% --------------------------- %

while xNow < xEnd
    
    prevSwitchIdx = memVarsIdx(someIdx);
    prevSwitchIdx = prevSwitchIdx(prevSwitchIdx > 0); % Removes any indices from forcing.
    
    Z = lagvals(xNow, prevSwitchIdx, lagsRVec, serialData); % Calculate memorisation variables from xNow.
    
    if INCL_FORCING
        
        Z2 = lagvals(xNow, [], lagsRVec, forcing);
        yNow = bdefun(Z, Z2);
        
    else
        
        yNow = bdefun(Z);
        
    end
    
    if ~isequal(yNow, sol.y(:, end)) % If a switch occurs at xNow (i.e. the values at xNow are different to our previously calculated point):
        
        sol.y = [sol.y, yNow]; % Append this to our solution.
        sol.x = [sol.x, xNow];
        
    end
    
    futPotSwitchIdx = potSwitches > xNow;
    xNow = min(potSwitches(futPotSwitchIdx)); % The next value of xNow.
    someIdx = find(potSwitches == xNow);
    
end

% --------------------------- %
% ------ End iteration ------ %
% --------------------------- %

% Post-processing of solution.

xNow = xEnd;

Z = lagvals(xNow, [], lagsRVec, serialData);

if INCL_FORCING
    
    Z2 = lagvals(xNow, [], lagsRVec, forcing);
    yNow = bdefun(Z, Z2);
    
else
    
    yNow = bdefun(Z);
    
end

sol.y = [sol.y, yNow];
sol.x = [sol.x, xNow];

sol.y = sol.y(:,lenHist + 1: end);
sol.x = sol.x(lenHist + 1: end);


% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

function Z = lagvals(tNow, prevSwitchIdx, lagsRVec, sol)
        
if tNow < sol.x(1) + max(lagsRVec)
    error('The value of tNow must be >= timeRVec(1) + max(lagsRVec).')
end
        
Z = false(size(sol.y, 1), length(lagsRVec)); % Initialise Z
        
for idxLag = 1:length(lagsRVec)
            
    idx = find(sol.x <= tNow - lagsRVec(idxLag));
            
    Z(:,idxLag) = sol.y(:, idx(end));
            
end
        
if ~isempty(prevSwitchIdx)
            
    idxMemDelays = memDelaysIdx(potSwitches == xNow);
    Z(:, idxMemDelays(idxMemDelays > 0)) = sol.y(:, prevSwitchIdx);
            
end
        
end

end