function [optNames, optPoss, optVals] = default_bdeopts
% function [optNames, optPoss, optVals] = default_bdeopts
% 
% DEFAULT_BDEOPTS   Set default options for the BDE solver.
%
% [optNames, optPoss, optVals] = default_bdeopts
%
% OUTPUTS
%
% optNames: List of options.
% optPoss: Type of corresponding input.
% optVals: Default values.
%
% INPUTS 
%
% None.
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
% Code review by Ozgur Akman, University of Exeter, 2019
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

optNames = {'MaxEvals'; 'MaxNumSwitches'}; % Option names
optPoss = {'positive scalar'; 'positive scalar'}; % Possible values
optVals = {Inf; Inf}; % Default values

end