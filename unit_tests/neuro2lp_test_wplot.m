function [solerr, solcell, solSercell] = neuro2lp_test_wplot
% function [solerr, solcell, solSercell] = neuro2lp_test_wplot
%
% NEURO2LP_TEST_WPLOT   Test function for the BDE solver - 2-loop Neurospora model of Akman et al. 
%                       This version also plots the parallel & serial solutions for comparison.
%
% [solerr, solcell, solSercell] = neuro2lp_test_wplot
% [~, ~, solSercell] = neuro2lp_test_wplot
%
% OUTPUTS 
%
% solerr: Vector of errors between the parallel and serial solutions across tests (calculated using continuous normalised Hamming distance).
%
% solcell/solSercell: Cells containing the solution structures returned by bdesolve/bdesolveserial, for each test (BOTH OPTIONAL). Each element
%                     of solcell/solSercell is a structure sol/solSer with the following common fields -
% sol(Ser).x: A vector containing the times of switch points.
% sol(Ser).y: A matrix with n rows where n is the number of state variables. Each column is the state following each switch.
% sol(Ser).history: The history used to generate the solution.
% 
% INPUTS
%
% None.
%
% MODEL
%
% The model equations are:
%
% frq(t) = (~FRQ1(t - tau_3) AND ~FRQ2(t - tau_4)) OR L(t - tau_5),
% FRQ1(t) = frq(t - tau_1),
% FRQ2(t) = frq(t - tau_2),
% L(t) = L(t - 24).
%
% Reference: Akman, O.E., Watterson, S., Parton, A., Binns, N., Millar, A.J. & Ghazal, P. J. Roy. Soc. Interface, 9(74) (2012).
%
% DEPENDENCIES 
%
% bdesolve.
% bdesolveserial.
% neuro2lp_bmod.
% bdecut.
% bdedist.
% bdeplot.
% bdemerge.
% 
% SEE ALSO 
%
% neuro2lp_test.
%
% -------------------------------------------------------------------------
%
% Written by Ozgur Akman, University of Exeter, 2019
% O.E.Akman@exeter.ac.uk
%

% Set the gates to 00111.

gvec = [0 0 1 1 1];

% Specify  model.

model = @(Z) neuro2lp_bmod(Z, gvec);

% Set the delays.

tau1 = 5.5; 
tau2 = 2;
tau3 = 6.5;
tau4 = 9;
tau5 = 9;

% Set the light period.

tau6 = 24;

lags = [tau1 tau2 tau3 tau4 tau5];
lagsall = [lags tau6];

% Set the range over which to integrate.

tEnd = 120;
tRange = [0 tEnd];

% Absorb the light period into the solution algorithms.

bdesolvetest = @(modlags, inthistory, intRange) bdesolve(model, [modlags tau6], inthistory, intRange);
bdesolveserialtest = @(modlags, solin, intRange) bdesolveserial(model, [modlags tau6], solin, intRange);

% ------------------------------------------------------------ %
% -------------------------- TEST 1 -------------------------- % 
% ------------------------------------------------------------ %
% Integrate from constant history in DD. 

% Set the (constant) history in DD.

history = [true; true; true; false];

% Solve the equations.

sol1 = bdesolvetest(lags, history, tRange);

% Solve the equations serially.

solSer1 = bdesolveserialtest(lags, sol1, tRange + max(lagsall));

% Extract the coincident regions of the two solutions.

[~, solCut1] = bdecut(sol1, solSer1.x(1));
[solSerCut1, ~] = bdecut(solSer1, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr1] = bdedist(solCut1, solSerCut1);

% Plot the parallel and serial solutions.

plotsolns(sol1, solSer1, tRange(end), 'DD from constant history');

% ------------------------------------------------------------ %
% -------------------------- TEST 2 -------------------------- % 
% ------------------------------------------------------------ %
% Integrate from limit cycle in DD. 

% Extract the last section of the oscillation to use as the initial history.

[~, sol1sec] = bdecut(sol1, sol1.x(end) - max(lagsall));
sol1sec.x = sol1sec.x - sol1sec.x(1);

% Solve the equations.

sol2 = bdesolvetest(lags, sol1sec, tRange(end));

% Solve the equations serially.

solSer2 = bdesolveserialtest(lags, sol2, [sol2.x(1) + max(lagsall) tRange(end) + max(lagsall)]);

% Extract the coincident regions of the two solutions.

[~, solCut2] = bdecut(sol2, solSer2.x(1));
[solSerCut2, ~] = bdecut(solSer2, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr2] = bdedist(solCut2, solSerCut2);

% Plot the parallel and serial solutions.

plotsolns(sol2, solSer2, tRange(end), 'DD from limit cycle');

% ------------------------------------------------------------ %
% -------------------------- TEST 3 -------------------------- % 
% ------------------------------------------------------------ %
% Integrate from DD limit cycle in 12:12 LD. 

% Extract the last section of the DD oscillation to use as the initial LD history.

[~, sol2sec] = bdecut(sol2, sol2.x(end) - max(lagsall));
sol2sec.x = sol2sec.x - sol2sec.x(1);

% Generate the forcing.

ldsol = ldcyc(sol2sec.x(end), 6, 18);

% Generate the history.

sol2sec.y(end,:) = [];
sol2secld = bdemerge(sol2sec, ldsol);

% Initial integration.

sol3 = bdesolvetest(lags, sol2secld, tRange(end));

% Extract the last section of the oscillation to use as the new LD history.

[~, sol3sec] = bdecut(sol3, sol3.x(end) - max(lagsall));
sol3sec.x = sol3sec.x - sol3sec.x(1);

% Solve the equations again.

sol3 = bdesolvetest(lags, sol3sec, tRange(end));

% Solve the equations serially.

solSer3 = bdesolveserialtest(lags, sol3, [sol3.x(1) + max(lagsall) tRange(end) + max(lagsall)]);

% Extract the coincident regions of the two solutions.

[~, solCut3] = bdecut(sol3, solSer3.x(1));
[solSerCut3, ~] = bdecut(solSer3, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr3] = bdedist(solCut3, solSerCut3);

% Plot the parallel and serial solutions.

plotsolns(sol3, solSer3, tRange(end), 'LD 12:12 from DD limit cycle');

% ------------------------------------------------------------ %
% -------------------------- TEST 4 -------------------------- % 
% ------------------------------------------------------------ %
% Integrate from 12:12 LD limit cycle in 6:18 LD. 

% Extract the last section of the 12:12 LD oscillation to use as the initial 6:18 LD history.

[~, sol3sec] = bdecut(sol3, sol3.x(end) - max(lagsall));
sol3sec.x = sol3sec.x - sol3sec.x(1);

% Generate the forcing.

ldsol = ldcyc(sol3sec.x(end), 8, 14);

% Generate the history.

sol3sec.y(end,:) = [];
sol3sec = bdemerge(sol3sec, ldsol);

% Initial integration.

sol4 = bdesolvetest(lags, sol3sec, tRange(end));

% Extract the last section of the oscillation to use as the new LD history.

[~, sol4sec] = bdecut(sol4, sol4.x(end) - max(lagsall));
sol4sec.x = sol4sec.x - sol4sec.x(1);

% Solve the equations again.

sol4 = bdesolvetest(lags, sol4sec, tRange(end));

% Solve the equations serially.

solSer4 = bdesolveserialtest(lags, sol4, [sol4.x(1) + max(lagsall) tRange(end) + max(lagsall)]);

% Extract the coincident regions of the two solutions.

[~, solCut4] = bdecut(sol4, solSer4.x(1));
[solSerCut4, ~] = bdecut(solSer4, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr4] = bdedist(solCut4, solSerCut4);

% Plot the parallel and serial solutions.

plotsolns(sol4, solSer4, tRange(end), 'LD 6:18 from LD 12:12 limit cycle');

% ------------------------------------------------------------ %
% -------------------------- TEST 5 -------------------------- % 
% ------------------------------------------------------------ %
% Integrate from 12:12 LD limit cycle in 16:8 LD. 

% Extract the last section of the 12:12 LD oscillation to use as the initial 6:18 LD history.

% Generate the forcing.

ldsol = ldcyc(sol3sec.x(end), 4, 20);

% Generate the history.

sol3sec.y(end,:) = [];
sol3sec = bdemerge(sol3sec, ldsol);

% Initial integration.

sol5 = bdesolvetest(lags, sol3sec, tRange(end));

% Extract the last section of the oscillation to use as the new LD history.

[~, sol5sec] = bdecut(sol5, sol5.x(end) - max(lagsall));
sol5sec.x = sol5sec.x - sol5sec.x(1);

% Solve the equations again.

sol5 = bdesolvetest(lags, sol5sec, tRange(end));

% Solve the equations serially.

solSer5 = bdesolveserialtest(lags, sol5, [sol5.x(1) + max(lagsall) tRange(end) + max(lagsall)]);

% Extract the coincident regions of the two solutions.

[~, solCut5] = bdecut(sol5, solSer5.x(1));
[solSerCut5, ~] = bdecut(solSer5, tRange(end));

% Calculate the distance between the solutions (i.e. the integration error).

[~, ~, solerr5] = bdedist(solCut5, solSerCut5);

% Plot the parallel and serial solutions.

plotsolns(sol5, solSer5, tRange(end), 'LD 18:6 from LD 12:12 limit cycle');

% Return the errors.

solerr = [solerr1 solerr2 solerr3 solerr4 solerr5];

% Gather together the parallel and serial solutions, if required.

if nargout > 1 
   
    solcell = {sol1 sol2 sol3 sol4 sol5};
    
    if nargout >2 
       
        solSercell = {solSer1 solSer2 solSer3 solSer4 solSer5};
        
    end
    
end

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction to generate a Boolean LD cycle.

function sol = ldcyc(tEnd, dawn, dusk)

ncycs = ceil(tEnd/24);
ftvec = [0 dawn dusk 24];    
fyvec = [0 1 0 0];
fyvec = logical(fyvec);
tvec = ftvec;
yvec = fyvec;
    
for k = 1:(ncycs - 1)        
    tvec=[tvec(1:end-1) ftvec + 24*k];
    yvec=[yvec(1:end-1) fyvec];        
end

in = find(tvec<tEnd);    
sol.x = [tvec(1:in(end)) tEnd];
sol.y = [yvec(1:in(end)) yvec(end)];
    
end

% Subfunction to plot the parallel and serial solutions.

function plotsolns(sol, solSer, tEnd, tstr)

figure; 
subplot(2, 1, 1);
bdeplot(sol);

if isfield(sol.history,'x')
    bdeplot(sol.history,[],[],'k');
end

subplot(2, 1, 2);
bdeplot(solSer);
bdeplot(solSer.history,[],[],'k');

for k=1:2
    subplot(2, 1, k);
    ax = axis;
    ax(1:2) = [0 tEnd];
    ax(3:4) = [-0.05 1.6];
    axis(ax);
    box on;
    xticks(0:4:tEnd); 
    yticks([0 0.1 0.2 0.3 1 1.1 1.2 1.3]); 
    yticklabels({'0', '0', '0', '0', '1', '1', '1', '1'}); 
    xlabel('Time, t (h)');
    ylabel('Expression (0/1)');
    set(gca, 'FontSize', 12);
    legend({'frq', 'FRQ1', 'FRQ2', 'L'}, 'Orientation', 'Horizontal', 'Location', 'North', 'FontSize', 18);
end

subplot(2, 1, 1);
title(['Parallel solution - ' tstr], 'FontSize', 20);
subplot(2, 1, 2);
title(['Serial solution - ' tstr], 'FontSize', 20);

end

end