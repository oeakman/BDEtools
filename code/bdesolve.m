function sol = bdesolve(bdefun, lags, history, xRange, varargin)
% function sol = bdesolve(bdefun, lags, history, xRange, varargin)
%
% BDESOLVE   Solve a system of Boolean Delay Equations.
%
% sol = BDESOLVE(bdefun, lags, history, xEnd)
% sol = BDESOLVE(bdefun, lags, history, xRange)
% sol = BDESOLVE(bdefun, lags, history, xEnd, forcing)
% sol = BDESOLVE(bdefun, lags, history, xEnd, forcing)
% sol = BDESOLVE(bdefun, lags, history, xEnd, forcing, options)
%
% OUTPUT 
%
% sol: A structure with the following fields -
% sol.x: A vector containing the times of switch points.
% sol.y: A matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% sol.solver: A string specifying the solver used.
% sol.bdefun: The system of equations that is solved.
% sol.lags: A vector of lag values.
% sol.numEvals: The number of times the equations were called.
% sol.history: The history used to generate the solution.
% sol.historyFlag: A logical value which specifies whether a switch occurs at the first point in the solution. This can indicate if the 
%                  history is incompatible with the solution. 
% sol.earlyTermFlag: A logical value which specifies whether the solution was terminated prematurely (e.g. if the maximum number of 
%                    switch points was exceeded).
% sol.forcing: This is only present if forcing was provided.
% 
% INPUTS
%
% bdefun: Handle to function that takes in a matrix of memorisation variables (each column corresponds to each value in <lags>) and returns 
%         a vector. If there is forcing, the function also takes in a matrix of memorisation variables for the forcing (see below).
%
% lags: A vector of lag values.
%
% history: Can be an n x 1 vector or a structure. If a structure, should contain at least two fields, history.x (vector containing times of 
%          switch points) and history.y (matrix where each row is the state after each switch). If a vector, this corresponds to a constant history.
% 
% xEnd/xRange: A scalar (if history is a structure) or a 1 x 2 vector (if history is a vector). If a scalar, this specifies the end point of integration 
%              and the start point is taken from the end point of the history. If a vector, it should contain the start and end points of integration.
%
% varargin{1} -
% forcing: A structure with two fields, forcing.x and forcing.y. forcing.x should span the same range of values as xRange. If forcing is 
% provided, then bdefun should take an extra argument, the memorisation variables for the forcing.
%
% varargin{2} -
% options: Structure that is the output of bdeopts, which can be used to specify extra options. 
%
% DEPENDENCIES
% 
% bdeopts.
%
% SEE ALSO 
%
% bdesolveserial (the serial solver).
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

% Set the final timepoint.

xEnd = xRange(end);

% Initialise output.

sol.bdefun= bdefun;
sol.history = history;
sol.numEvals = 0;
sol.solver = 'bdesolve';
sol.lags = lags;
sol.earlyTermFlag = false;

% Initialise history.

if ~isstruct(history)
    
    if size(history, 2) > 1
        
        error('History must be specified as a structure or a column vector.')
        
    end
    
    history1 = struct;
    history1.x = [xRange(1) - max(lags), xRange(1)];
    history1.y = [history history];
    history = history1;
    
end

if numel(xRange) == 1 && (islogical(history) || isnumeric(history))
    
    error('If history is a constant vector then an initial x must be provided, i.e. a two-element vector, xRange, is required.')
    
end

if isstruct(history)
    
    if ~isequal(numel(history.x), numel(history.y(1,:)))
        error('The number of elements in history.x should match the number of columns in history.y')
    end
    
    if ~islogical(history.y)
        
        if any(history.y ~= 0 || history.y ~= 1)
            
            error('History is not logical. It is recommended to provide the history as a logical vector.')
            
        end
        
        history.y = logical(history.y);
        
    end
    
else
    
    if ~islogical(history)
        
        if any(history ~= 0 || history ~= 1)
            
            error('History is not logical. It is recommended to provide the history as a logical vector.')
            
        end
        
        history = logical(history);
        
    end
    
    if size(history, 2) > 1
        
        error('History must be specified as a structure or a column vector.')
        
    end
    
    history1 = struct;
    history1.x = [xRange(1) - max(lags), xRange(1)];
    history1.y = [history history];
    history = history1;
    
end

% Initialise forcing.

if nargin > 4
    
    if ~isempty(varargin{1})
        
        INCL_FORCING = true;
        forcing = varargin{1};
        
        if ~isstruct(forcing)
            
            error('Forcing must be specified as a structure.')
            
        end
        
        if ~islogical(forcing.y)
            
            if ~all(forcing.y == 0 | forcing.y == 1)
                
                error('Forcing.y is not logical. It is recommended to provide forcing.y as a logical vector.')
                
            end
            
            forcing.y = logical(forcing.y);
            
        end
        
        if (forcing.x(1) ~= history.x(1)) || (forcing.x(end) ~= xEnd)
            
            error('The forcing should span the same range as the model prediction and history.')
            
        end
        
        sol.forcing = forcing;
        
        [numForceVars, numForceSwitches] = size(forcing.y);
        
    else
        
        INCL_FORCING = false;
        
    end
    
else
    
    INCL_FORCING = false;
    
end

% Set BDE solver options.

if nargin > 5 % If an options structure is provided.
    
    my_options = varargin{2};
    
    options = bdeopts; % Calling bdeopts with no arguments returns the default options.
    
    optNames = fieldnames(options);
    
    idxOpts = structfun(@(x) ~isempty(x), my_options); % Determine which options have been specified by the user.
    
    for k=1:length(idxOpts)    
        options.(optNames{idxOpts(k)}) = my_options.(optNames{idxOpts(k)}); % Set the values of these variables to the user-defined values.    
    end
        
else
    
    options = bdeopts;
    
end

sol.bdeoptions=options;

numLags = numel(lags);

numVars = numel(history.y(:,1));
lenHist = numel(history.x);

candBatchSize = 10000;  % Initialise this many entries for the vector containing candidate switches and others, to
                        % avoid dynamic memory allocation.
solBatchSize = 1000;    % Initialise this many columns for the solution.
solBatchEnd = solBatchSize;

candidates = nan(1, candBatchSize); % The values of candidate switches.
candSwitchIdx = zeros(1, candBatchSize, 'uint64'); % The index of which switch in sol.x that generated each candidate.
candDelaysIdx = zeros(1, candBatchSize, 'uint64'); % The index of which delay that generated each candidate.
candDelaysSum = zeros(numLags + 1, candBatchSize, 'uint64'); % First row is an index into the point in the history.
candVarsIdx = false(numVars, candBatchSize); % Record which variable has changed value to produce each candidate.
memVars = false(numVars, numLags);
solDelaysSum = zeros(numLags + 1, solBatchSize, 'uint64'); % A matrix where the first row indexes the point in
                                                           % the history where each switch "originates" and the remaining rows give
                                                           % the number of each delay that can be added to the original history point
                                                           % to get the current point.

if INCL_FORCING
    
    solDelaysSumForce = zeros(numLags + 1, solBatchSize, 'uint64'); % Like solDelaysSum but the first row indexes
                                                                    % the switches in the forcing vector that result in a candidate.
    candDelaysSumForce = zeros(numLags + 1, candBatchSize, 'uint64');
    candVarsForceIdx = false(numForceVars, candBatchSize);
    candSwitchForceIdx = zeros(1, candBatchSize, 'uint64');
    xNowIdx = 1;
    
end

candCount = 1; % Counter for the vector of candidates.

if INCL_FORCING
    
    xNowForceIdx = 2;
    
    % Calculate candidate switch points caused by the forcing.
    
    for idxForceSwitches = 2:numForceSwitches
        
        solDelaysSumForce(1, xNowForceIdx) = idxForceSwitches;
        
        [candDelaysIdx, candDelaysSumForce, candVarsIdx, candCount, candidates, candSwitchForceIdx] = ...
            calc_future_switches(forcing, candidates, candDelaysIdx, candDelaysSumForce, candVarsIdx, solDelaysSumForce, candCount, numLags, numVars, [], xNowForceIdx, candBatchSize, lags, candSwitchForceIdx);
        
        xNowForceIdx = xNowForceIdx + 1;
        
        change = xor(forcing.y(:, idxForceSwitches - 1), forcing.y(:, idxForceSwitches));
        candVarsForceIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
        
        candCount = candCount + numLags;
        
    end
    
end

xNowIdx = 1;

for idx_hist = 1:lenHist - 2 % Initially, add the relevant lags to each of the switch points in the history.
    
    xNow = history.x(idx_hist + 1);
    
    xNowIdx = xNowIdx + 1;
    
    solDelaysSum(1, xNowIdx) = idx_hist + 1;
    
    [candDelaysIdx, candDelaysSum, candVarsIdx, candCount, candidates, candSwitchIdx] = ...
        calc_future_switches(history, candidates, candDelaysIdx, candDelaysSum, candVarsIdx, solDelaysSum, candCount, numLags, numVars, xNow, xNowIdx, candBatchSize, lags, candSwitchIdx);
    
    change = xor(history.y(:, idx_hist), history.y(:, idx_hist + 1));
    candVarsIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
    
    candCount = candCount + numLags;
    
end


xNow = history.x(end);
for idxLag1 = 1:numLags
    
    idxH = find(history.x <= xNow - lags(idxLag1));
    
    memVars(:, idxLag1) = history.y(:, idxH(end));
    
end

sol.x = [history.x(1: end - 1), nan(1, solBatchSize)];
sol.y = [history.y(:, 1: end - 1), false(numVars, solBatchSize)];

% Need to check the last point in the history as a switch may be necessary here (to patch the history with the prediction):
xNow = history.x(end);
xNowIdx = lenHist;

Z = lagvals(xNow, xNowIdx, [], lags, sol, numVars, numLags, candDelaysIdx, candidates);
memVars = Z;

if INCL_FORCING
    
    Z2 = lagvals(xNow, xNowIdx, [], lags, forcing, numForceVars, numLags, candDelaysIdx, candidates);
    yNow = bdefun(Z, Z2);
    
    memVarsF = Z2;
    
else
    
    yNow = bdefun(Z);
    
end

sol.numEvals = sol.numEvals + 1;

sol.x = [history.x(1: end - 1), xNow, nan(1, solBatchSize)]; % First point in solution is last point in history.
sol.y = [history.y(:, 1: end - 1), yNow, false(numVars, solBatchSize)]; 


if ~isequal(history.y(:, end-1), yNow)
    
    % If a switch occurs at xNow (the last point in the history):
    
    solDelaysSum(1, xNowIdx) = lenHist;
    
    % Calculate future potential switch points as a result of this switch.
    
    if solDelaysSum(1, xNowIdx) ~= 0
        
        [candDelaysIdx, candDelaysSum, candVarsIdx, candCount, candidates, candSwitchIdx] = ...
            calc_future_switches(history, candidates, candDelaysIdx, candDelaysSum, candVarsIdx, solDelaysSum, candCount, numLags, numVars, xNow, xNowIdx, candBatchSize, lags, candSwitchIdx);
        
        change = xor(history.y(:, end-1), yNow);
        candVarsIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
        
        candCount = candCount + numLags;
        
    else
        
        [candDelaysIdx, candDelaysSumForce, candVarsIdx, candCount, candidates, candSwitchIdx] = ...
            calc_future_switches(forcing, candidates, candDelaysIdx, candDelaysSumForce, candVarsIdx, solDelaysSumForce, candCount, numLags, numVars, xNow, xNowIdx, candBatchSize, lags, candSwitchIdx);
        
        change = xor(history.y(:, end-1), yNow);
        candVarsIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
        
        candCount = candCount + numLags;
        
    end
    
    if ~isequal(yNow, history.y(:, end))
        % If the switch is one that must be imposed (i.e. not equal to the last point in the history provided by the user):
        
        sol.historyFlag = true; % historyFlag is assigned the value 1 if it is necessary
                                % to incorporate a switch to patch the end of the history to the start
                                % of the prediction.        
    else
        
        sol.historyFlag = false;
        
    end
    
else
    
    sol.historyFlag = false;
    
end

futureCandidateIdx = candidates > xNow;

if ~any(futureCandidateIdx > 0)
    
    xNow = xEnd;
    
else
    
    xNow = min(candidates(futureCandidateIdx)); % The next value of xNow.
    currCandIdx = find(candidates == xNow);
    
    if isempty(currCandIdx)
        
        error('Something has gone wrong!') % Sanity check.
        
    end
    
end

term = false; % term = true if the solution is being terminated.

% --------------------------- %
% ----- Begin iteration ----- %
% --------------------------- %

while xNow < xEnd 
    
    for idxCand = 1:numel(currCandIdx)
        
        if ~ all(candVarsIdx(:, currCandIdx(idxCand)) == 0) % This should only ever be not satisfied when including forcing for the case where the memorisation variables for
                                                            % forcing switch but the other memorisation variables do not.
            
            memVars(candVarsIdx(:, currCandIdx(idxCand)), candDelaysIdx(currCandIdx(idxCand))) = ~memVars(candVarsIdx(:, currCandIdx(idxCand)), candDelaysIdx(currCandIdx(idxCand)));
            
        end
        
        if INCL_FORCING && ~ all(candVarsForceIdx(:, currCandIdx(idxCand)) == 0) % If there is a switch in the memorisation variables for forcing.
            
            memVarsF(candVarsForceIdx(:, currCandIdx(idxCand)), candDelaysIdx(currCandIdx(idxCand))) = ~memVarsF(candVarsForceIdx(:, currCandIdx(idxCand)), candDelaysIdx(currCandIdx(idxCand)));
            
        end
        
    end
    
    prevSwitchIdx = candSwitchIdx(currCandIdx);
    prevSwitchIdx = prevSwitchIdx(prevSwitchIdx > 0); % Removes any indices from forcing.
    
    Z = memVars;
    
    if INCL_FORCING
        
        Z2 = memVarsF;
        yNow = bdefun(Z, Z2);
        
    else
        
        yNow = bdefun(Z);
        
    end
    
    sol.numEvals = sol.numEvals + 1;
    
    if sol.numEvals >= options.MaxEvals || xNowIdx >= lenHist + options.MaxNumSwitches
        
        term = true;
        
        sol.earlyTermFlag = true;
        
        xNowIdx = xNowIdx + 1;
        
        if sol.numEvals >= options.MaxEvals
            
            disp_message = 'Terminating solution: Maximum number of evaluations reached.';
            
        elseif xNowIdx >= options.MaxNumSwitches
            
            disp_message = 'Terminating solution: Maximum number of switches reached.';
            
        end
        
        break;
        
    end    
    
    if ~isequal(yNow, sol.y(:, xNowIdx)) % If a switch occurs at xNow (i.e. the values at xNow are different to our previously calculated point):
        
                
        if xNowIdx > solBatchEnd
            
            [sol, solBatchEnd] = update_sol_storage(sol, numVars, solBatchSize, solBatchEnd);
            
        end
        
        xNowIdx = xNowIdx + 1;
        
        sol.y(:, xNowIdx) = yNow; % Append this to our solution.
        sol.x(xNowIdx) = xNow;
        for cIdx = 1:numel(currCandIdx)
            if candSwitchIdx(currCandIdx(cIdx)) ~= 0 && solDelaysSum(1, candSwitchIdx(currCandIdx(cIdx))) ~= 0
                
                solDelaysSum(:, xNowIdx) = [solDelaysSum(1, candSwitchIdx(currCandIdx(cIdx))); ...
                    solDelaysSum(2:end, candSwitchIdx(currCandIdx(cIdx))) + uint64((1:numLags == candDelaysIdx(currCandIdx(cIdx)))')];
                % NOTE: Indexing the first element of currCandIdx is necessary when there is more than one way to sum up to the current point.
                
            else
                
                solDelaysSumForce(:, xNowForceIdx) = [solDelaysSumForce(1, candSwitchForceIdx(currCandIdx(cIdx))); ...
                    solDelaysSumForce(2:end, candSwitchForceIdx(currCandIdx(cIdx))) + uint64((1:numLags == candDelaysIdx(currCandIdx(cIdx)))')];
                
            end
        end
        
        if solDelaysSum(1, xNowIdx) ~= 0
            
            [candDelaysIdx, candDelaysSum, candVarsIdx, candCount, candidates, candSwitchIdx] = ...
                calc_future_switches(history, candidates, candDelaysIdx, candDelaysSum, candVarsIdx, solDelaysSum, candCount, numLags, numVars, xNow, xNowIdx, candBatchSize, lags, candSwitchIdx);
            
            change = xor(sol.y(:, xNowIdx - 1), sol.y(:, xNowIdx));
            candVarsIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
            
            candCount = candCount + numLags;
            
        else
            
            [candDelaysIdx, candDelaysSumForce, candVarsIdx, candCount, candidates, candSwitchForceIdx] = ...
                calc_future_switches(forcing, candidates, candDelaysIdx, candDelaysSumForce, candVarsIdx, solDelaysSumForce, candCount, numLags, numVars, xNow, xNowForceIdx, candBatchSize, lags, candSwitchForceIdx);
            
            change = xor(sol.y(:, xNowIdx - 1), sol.y(:, xNowIdx));
            candVarsIdx(:, candCount:candCount + numLags - 1) = repmat(change, 1, numLags);
            
            candCount = candCount + numLags;
            
        end
                        
    end
    
    futureCandidateIdx = candidates > xNow;
    
    if sum(futureCandidateIdx) == 0
        
        xNow = xEnd;
        
    else
        
        xNow = min(candidates(futureCandidateIdx)); % The next value of xNow.
        
        currCandIdx = find(candidates == xNow); 
        
    end
    
end 

% --------------------------- %
% ------ End iteration ------ %
% --------------------------- %

% Post-processing of solution.

if term
    
    if xNowIdx > solBatchEnd
        
        [sol, solBatchEnd] = update_sol_storage(sol, numVars, solBatchSize, solBatchEnd);
        
    end
    
    % Include the last evaluated point:
    
    sol.y(:, xNowIdx) = yNow;
    sol.x(xNowIdx) = xNow;
    
    disp(disp_message);
    
else
    
    % Include a point for xEnd:
    
    xNow = xEnd;
    
    Z = lagvals(xNow, xNowIdx, [], lags, sol, numVars, numLags, candDelaysIdx, candidates); % Calculate memorisation variables from tNow.
    
    if INCL_FORCING
        
        Z2 = lagvals(xNow, numel(forcing.x), [], lags, forcing, numForceVars, numLags, candDelaysIdx, candidates);
        yNow = bdefun(Z, Z2);
        
    else
        
        yNow = bdefun(Z);
        
    end
    
    sol.numEvals = sol.numEvals + 1;
    
    xNowIdx = xNowIdx + 1;
    
    sol.y(:, xNowIdx) = yNow;
    sol.x(xNowIdx) = xNow;
    
end

sol.x = sol.x(1:xNowIdx);
sol.y = sol.y(:, 1:xNowIdx);

sol.y = sol.y(:,lenHist: end);
sol.x = sol.x(lenHist: end);


% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction to calculate future switches.

function [candDelaysIdx, candDelaysSum, candVarsIdx, candCount, candidates, candSwitchIdx] = ...
                calc_future_switches(history, candidates, candDelaysIdx, candDelaysSum, candVarsIdx, solDelaysSum, candCount, numLags, numVars, xNow, xNowIdx, candBatchSize, lags, candSwitchIdx)
                
candDelaysIdx(candCount:candCount + numLags - 1) = 1:numLags;
candDelaysSum(:, candCount:candCount + numLags - 1) = solDelaysSum(:, xNowIdx) + [zeros(1, numLags, 'uint64'); eye(numLags, 'uint64')];
                
newTimes = history.x(candDelaysSum(1, candCount:candCount + numLags - 1)) + sum(double(candDelaysSum(2:end, candCount:candCount + numLags - 1)) .* repmat(lags', 1, numLags), 1);
        
if candCount > candBatchSize && candCount <= candBatchSize + numLags
            
[candidates, candDelaysIdx, candSwitchIdx, candDelaysSum, candVarsIdx] = ...
                update_candidate_storage(candidates, candDelaysIdx, candSwitchIdx, candDelaysSum, candVarsIdx, candBatchSize, numLags, numVars);
            
end
        
candidates(candCount:candCount + numLags - 1) = newTimes(:)';
candSwitchIdx(candCount:candCount + numLags - 1) = xNowIdx * ones(1, numel(newTimes));
               
end

% --------------------------------------- %

% Subfunction to generate matrix of lag values.

function Z = lagvals(xNow, xNowIdx, prevSwitchIdx, lags, sol, numVars, numLags, candDelaysIdx, candidates)
        
Z = false(numVars, numLags); % Initialise Z.
        
for idxLag = 1:numLags
            
    idx = find(sol.x(1:xNowIdx) <= xNow - lags(idxLag));
            
    Z(:, idxLag) = sol.y(:, idx(end));
            
end
        
if ~isempty(prevSwitchIdx)
            
    idxMemDelays = candDelaysIdx(candidates == xNow);
    Z(:, idxMemDelays(idxMemDelays > 0)) = sol.y(:, prevSwitchIdx);
            
end
        
end

% --------------------------------------- %

% Subfunction used by calc_future_switches.

function [candidates, candDelaysIdx, candSwitchIdx, candDelaysSum, candVarsIdx] = ...
            update_candidate_storage(candidates, candDelaysIdx, candSwitchIdx, candDelaysSum, candVarsIdx, candBatchSize, numLags, numVars)
        
% Add extra columns to matrices if we have reached our limit. We could remove many candidates at this stage since they are not
% required but for simplicity, candidates have been kept for the moment.
        
candidates = [candidates, nan(1, candBatchSize)];
candDelaysIdx = [candDelaysIdx, zeros(1, candBatchSize, 'uint64')];
candSwitchIdx = [candSwitchIdx, zeros(1, candBatchSize, 'uint64')];
candDelaysSum = [candDelaysSum, zeros(numLags + 1, candBatchSize, 'uint64')];
candVarsIdx = [candVarsIdx, false(numVars, candBatchSize)];
        
end

% --------------------------------------- %

% Subfunction to update output structure.

function [sol, solBatchEnd] = update_sol_storage(sol, numVars, solBatchSize, solBatchEnd)
        
sol.x = [sol.x, nan(1, solBatchSize)];
sol.y = [sol.y, false(numVars, solBatchSize)];
solBatchEnd = solBatchEnd + solBatchSize;
        
end

end