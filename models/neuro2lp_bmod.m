function yNow = neuro2lp_bmod(Z, g)
% function yNow = neuro2lp_bmod(Z, g)
%
% NEURO2LP_BMOD   BDE version of 2-loop Neurospora model from Akman et al. 
%
% OUTPUT
%
% yNow: Column vector of the current values of the model variables. 
%
% INPUTS
%
% Z: A 4x6 matrix of memorisation variables.
% g: A Boolean vector determining whether a particular single-input logic gate returns its input or performs the
%    NOT operation (i.e. whether a reaction is activating or inhibitory), and whether a double-input logic gate peforms 
%    an AND or an OR operation.
%
% MODEL
%
% The model equations are:
%
% frq(t)= G2(G1(FRQ1(t - tau_3), g_3), G1(FRQ2(t - tau_4), g_4), g_5) OR L(t - tau_5),
% FRQ1(t) = G1(frq(t - tau_1), g_1),
% FRQ2(t) = G1(frq(t - tau_2), g_2),
% L(t) = L(t - tau_6),
%
% where G1(X, 0) = X, G1(X, 1) = NOT X, G2(X, Y, 0) = X OR Y, G2(X, Y, 1) = X AND Y and tau_6 is the period of the light input (normally, tau_6 = 24).
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

mRNA1 = applygate1(Z(2, 3), g(3));
mRNA2 = applygate1(Z(3, 4), g(4));

yNow = [applygate2(mRNA1, mRNA2, g(5)) | Z(4, 5);
        applygate1(Z(1, 1), g(1));
        applygate1(Z(1, 2), g(2));
        Z(4, 6)];

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

% Subfunction defining G2.

function bitout = applygate2(bitin1, bitin2, inputgate)

if inputgate
    bitout = bitin1 & bitin2; % AND gate.
else
    bitout = bitin1 | bitin2; % OR gate.    
end

end

end