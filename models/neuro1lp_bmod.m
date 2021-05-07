function yNow = neuro1lp_bmod(Z, g)
% function yNow = neuro1lp_bmod(Z, g)
%
% NEURO1LP_BMOD   BDE version of 1-loop Neurospora model from Akman et al. 
%
% OUTPUT
%
% yNow: Column vector of the current values of the model variables. 
%
% INPUTS
%
% Z: A 3x4 matrix of memorisation variables.
% g: A Boolean vector determining whether a particular logic gate returns its input or performs the
%    NOT operation (i.e. whether a reaction is activating or inhibitory).
%
% MODEL
%
% The model equations are:
%
% frq(t) = G1(FRQ(t - tau_2), g_2) OR L(t - tau_3),
% FRQ(t) = G1(frq(t - tau_1), g_1),
% L(t) = L(t - tau_4),
%
% where G1(X, 0) = X and G1(X, 1) = NOT X and tau_4 is the period of the light input (normally, tau_4 = 24). 
%
% Reference: Akman, O.E., Watterson, S., Parton, A., Binns, N., Millar, A.J. & Ghazal, P. J. Roy. Soc. Interface, 9(74) (2012).
%
% DEPENDENCIES 
%
% None.
% 
% -------------------------------------------------------------------------
%
% Written by Ozgur Akman & Kevin Doherty, University of Exeter, 2017
% O.E.Akman@exeter.ac.uk
% K.Doherty@exeter.ac.uk
%
% Code review and edits by Ozgur Akman, 2019
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

yNow = [applygate1(Z(2, 2), g(2)) | Z(3, 3);
        applygate1(Z(1, 1), g(1));
        Z(3, 4)];

% --------------------------------------- %
% ----------- SUBFUNCTIONS--------------- %
% --------------------------------------- %

% Subfunction defining G1.

function bitout = applygate1(bitin, inputgate)

if inputgate 
    bitout = ~bitin; % NOT gate.
else
    bitout = bitin;  % Identity gate.
end

end

end