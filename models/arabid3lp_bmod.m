function yNow = arabid3lp_bmod(Z, g)
% function yNow = arabid3lp_bmod(Z, g)
%
% ARABID3LP_BMOD   BDE version of 3-loop Arabidopsis model from Akman et al. 
%
% OUTPUT
%
% yNow: Column vector of the current values of the model variables. 
%
% INPUTS
%
% Z: A 9x13 matrix of memorisation variables.
% g: A Boolean vector determining whether a particular single-input logic gate returns its input or performs the
%    NOT operation (i.e. whether a reaction is activating or inhibitory), and whether a double-input logic gate peforms 
%    an AND or an OR operation.
%
% MODEL
%
% The model equations are:
%
% LHY(t) = G2(G1(X(t - tau_3), g_3), G1(PRR(t - tau_8), g_10), g_11) AND L1(t - tau_9),
% TOC1(t) = G2(G1(LHY(t - tau_1), g_1), G1(Y(t - tau_6), g_6), g_8),
% X(t) = G1(TOC1(t - tau_2), g_2),
% Y(t) = G2(G1(LHY(t - tau_4), g_4), G1(TOC1(t - tau_5), g_5), g_7) AND (L2(t - tau_10) OR L3(t - tau_11)),
% PRRR(t) = G1(LHY(t - tau_7), g_9) AND L4(t - tau_12),
% L1(t) = L1(t - tau_13),
% L2(t) = L2(t - tau_13),
% L3(t) = L3(t - tau_13),
% L4(t) = L4(t - tau_13),
%
% where G1(X, 0) = X, G1(X, 1) = NOT X, G2(X, Y, 0) = X OR Y, G2(X, Y, 1) = X AND Y and tau_13 is the period of the light input (normally, tau_13 = 24).
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

yNow = [applygate2(applygate1(Z(3, 3), g(3)), applygate1(Z(5, 8), g(10)), g(11)) & Z(6, 9);
        applygate2(applygate1(Z(1, 1), g(1)), applygate1(Z(4, 6), g(6)), g(8));
        applygate1(Z(2, 2), g(2));
        applygate2(applygate1(Z(1, 4), g(4)), applygate1(Z(2, 5), g(5)), g(7)) & (Z(7, 10) | Z(8, 11));
        applygate1(Z(1, 7), g(9)) & Z(9, 12);        
        Z(6, 13);
        Z(7, 13);
        Z(8, 13);
        Z(9, 13)];

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