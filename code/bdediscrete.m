function sol = bdediscrete(x, y, T, varargin)
% function sol = bdediscrete(x, y, T, varargin)
%
% BDEDISCRETE   Discretise a timeseries by thresholding.
%
% sol = bdediscrete(x, y, T)
% sol = bdediscrete(x, y, T, interpfun)
%
% OUTPUT
%
% sol: A structure with the following fields - 
% sol.x: A vector with the times of switch points.
% sol.y: A Boolean matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% 
% INPUTS
%
% x: Vector of time points.
% 
% y: Matrix of expression levels. Has n rows where each row is the timeseries associated with a different variable. 
% 
% T: Vector of thresholds, given as proportions of the ranges of the expression data. Each element should be between 0 and 1 where 0 corresponds
%    to the threshold being equal to the minimum expression level and 1 corresponds to the threshold being equal to the maximum expression level.
% 
% varargin{1} -
% interpfun: A custom interpolation function. The default is linear interpolation between successive points. 
%
% DEPENDENCIES 
%
% bdereduce.
%
% -------------------------------------------------------------------------
%
% Written by Ozgur Akman & Kevin Doherty, University of Exeter, 2017
% O.E.Akman@exeter.ac.uk
% k.doherty@exeter.ac.uk
%
% Code review by Ozgur Akman, University of Exeter, 2019
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

if nargin > 3    
    interpfun = varargin{1};    
else     
    interpfun = @lin_interp;    
end

Ttrans = T';

yNorm = (y - min(y, [], 2))./(max(y, [], 2) - min(y, [], 2)); % Normalise the time series between 0 and 1.

yThresh = yNorm >= Ttrans; % Boolean matrix with values determined by whether expression level is above (1) or below (0) threshold.

yThresh([false(size(y, 1), 1), yNorm(:, 2:end-1) == Ttrans & (yNorm(:, 1:end-2) < Ttrans == yNorm(:, 3:end) < Ttrans), false(size(y, 1), 1)]) = false; 
% Ignore infinitessimally small switches. This line works by removing true values from yThresh if the value of a data point is exactly equal to a threshold 
% at a point but the data does not cross the threshold.

switchesMat = diff(yThresh, [], 2); % Find where there are switches from 0 to 1 (= 1) and from 1 to 0 (= -1).

[~, switchIdxCVec] = find(switchesMat ~= 0); % Get the indices for where the switches occur (an index is repeated if it occurs in several variables).

switchIdxCVec = switchIdxCVec(:); % This condition is included to account for the fact that the output of find can be a row or a column vector.

switchIdxUniqueCVec = unique(switchIdxCVec); % Remove duplicates.

% Interpolate across switch points.

y1 = yNorm(:, switchIdxUniqueCVec); 
y2 = yNorm(:, switchIdxUniqueCVec + 1);
x1 = x(:, switchIdxUniqueCVec);
x2 = x(:, switchIdxUniqueCVec + 1);

[xI, yI] = interpfun(x1, y1, x2, y2); % Interpolate.

sol.x = xI;
sol.y = yI;


sol = bdereduce(sol); % Remove any duplicate switches.

    function [xI, yI] = lin_interp(x1, y1, x2, y2)
        
        % Note: accepts normalised y values. 
        
        xInterpMat = x1 + (x2 - x1) ./ (y2 - y1) .* (Ttrans - y1); % Calculates the time points where the linear interpolation between the points (x1,y1) 
                                                                   % and (x2,y2) is equal to the threshold value. Note: This also calculates values where 
                                                                   % the threshold does not lie between y1 and y2.
        
        idxPointsToKeep = (xInterpMat >= x1 & xInterpMat <= x2); % Indices of the points where the interpolant crosses the threshold between 
                                                                 % consecutive y values. We wish to remove all other points.
        
        
        % We still need to remove infinitessimal switches (although some may have been removed in the parent function, they may be reintroduced here if 
        % they occur around the same time as another switch due to the fact that we vectorised the calculation of xInterpMat):
        
        switchAtBoundaryIdx1 = xInterpMat == x1;
        switchAtBoundaryIdx2 = xInterpMat == x2;
        [~, b1] = find(switchAtBoundaryIdx1);
        [~, b2] = find(switchAtBoundaryIdx2);
        boundSwitchIdx1 = switchIdxUniqueCVec(any(switchAtBoundaryIdx1, 1));
        boundSwitchIdx2 = switchIdxUniqueCVec(any(switchAtBoundaryIdx2, 1)); 

        if ~isempty(boundSwitchIdx1) 
            idxBla = boundSwitchIdx1 > 1;
            boundSwitchIdx1 = boundSwitchIdx1(idxBla); % Leave infinitessimal switches if they occur in the first or last point of the data.
            b1_2 = b1(idxBla);
            switchAtBoundaryIdx11 = switchAtBoundaryIdx1(:, b1_2);
            if ~isempty(boundSwitchIdx1)
                ignoreIdx1 = yThresh(:, boundSwitchIdx1 - 1) == yThresh(:, boundSwitchIdx1 + 1); % Indices of infinitessimal switches.
                
                for i = 1:size(switchAtBoundaryIdx11, 2) 
                    idxPointsToKeep(logical(ignoreIdx1(:, i) .* switchAtBoundaryIdx11(:, i)), b1_2(i)) = false;
                end
            end
        end
        
        if ~isempty(boundSwitchIdx2) 
            idxBla = boundSwitchIdx2 < numel(x) - 1;
            boundSwitchIdx2 = boundSwitchIdx2(idxBla); % Leave infinitessimal switches if they occur in the first or last point of the data.
            b2_2 = b2(idxBla);
            switchAtBoundaryIdx12 = switchAtBoundaryIdx2(:, b2_2);
            if ~isempty(boundSwitchIdx2)
                ignoreIdx2 = yThresh(:, boundSwitchIdx2) == yThresh(:, boundSwitchIdx2 + 2); % Indices of infinitessimal switches.
                
                for i = 1:size(switchAtBoundaryIdx12, 2) 
                    idxPointsToKeep(logical(ignoreIdx2(:, i) .* switchAtBoundaryIdx12(:, i)), b2_2(i)) = false;
                end
            end
        end
        
        numSwitches = sum(idxPointsToKeep, 1); % Calculate the number of switches between each set of time points.
        
        xInterpGoodMat = xInterpMat;
        xInterpGoodMat(~idxPointsToKeep) = NaN; % Set any unwanted values to NaN.
        
        xInterpDesiredRVec = sort(xInterpMat(idxPointsToKeep)); % Only keep the desired calculated values and order them.
        xInterpDesiredRVec = xInterpDesiredRVec(:)'; % Ensure vector is a row vector.
        
        if ~isempty(xInterpDesiredRVec) && isequal(x(1), xInterpDesiredRVec(1))
            
            xI = [xInterpDesiredRVec, x(end)];
            
        else
            
            xI = [x(1), xInterpDesiredRVec, x(end)]; % Combine all desired switch times (only those where x1 < xInterpMat < x2) in one vector 
                                                     % and include initial and end points.
        
        end
        
        yI = yThresh(:, [1; switchIdxCVec + 1; end]); % The y values after each switch.
        
        % At this point yI contains values for each switch, but some of these can occur at the same time. We now recalculate where there
        % is more than one switch at a time point:
        
        numChanges1 = sum(~isnan(xInterpGoodMat), 1); % The number of variables that change value at each switch point.
        numChanges2 = sum(abs(diff(yI, [], 2)), 1);
        
        idxMultChanges1 = find(numChanges1 > 1); % Indices where there is more than one change.
        idxMultChanges2 = find(numChanges2 > 1);
        [numRows, ~] = size(yThresh);
        
        if ~isempty(idxMultChanges1) && (numel(idxMultChanges1) == numel(idxMultChanges2))
            for i = 1:numel(idxMultChanges1)
                idxSort = sortrows([xInterpGoodMat(:, idxMultChanges1(i)), (1:numel(xInterpGoodMat(:, 1)))']);
                for j = 1: numSwitches(idxMultChanges1(i))
                    idx2Change = false(numRows, 1); 
                    idx2Change(idxSort(j, 2)) = true;
                    yI(:, idxMultChanges2(i) + j) = xor(yI(:, idxMultChanges2(i) + j - 1), idx2Change);
                end
            end
        end
        
        % Remove entries where there is a switch in more than one variable at the same time.
        
        yI(:, diff(xI) == 0) = [];
        xI = unique(xI);
        
    end

end