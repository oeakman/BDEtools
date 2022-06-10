function yNow = arabid2lpnoX_bmod(Z, g)
% function yNow = arabid2lpnoX_bmod(Z, g)
%
% ARABID2LPNOX_BMOD   BDE version of 2-loop Arabidopsis model from Akman et al, with gene X removed.
%
% OUTPUT
%
% yNow: Column vector of the current values of the model variables. 
%
% INPUTS
%
% Z: A 6x9 matrix of memorisation variables.
% g: A Boolean vector determining whether a particular single-input logic gate returns its input or performs the
%    NOT operation (i.e. whether a reaction is activating or inhibitory), and whether a double-input logic gate peforms 
%    an AND or an OR operation.
%
% MODEL
%
% The model equations are:
%
% LHY(t) = G1(TOC1(t - tau_2), g_2) AND L1(t - tau_6),
% TOC1(t) = G2(G1(LHY(t - tau_1), g_1), G1(Y(t - tau_5), g_5), g_7),
% Y(t) = G2(G1(LHY(t - tau_3), g_3), G1(TOC1(t - tau_4), g_4), g_6) AND (L2(t - tau_7) OR L3(t - tau_8)),
% L1(t) = L1(t - tau_9),
% L2(t) = L2(t - tau_9),
% L3(t) = L2(t - tau_9),
%
% where G1(X, 0) = X, G1(X, 1) = NOT X, G2(X, Y, 0) = X OR Y, G2(X, Y, 1) = X AND Y tau_9 is the period of the light input (normally, tau_9 = 24).
%
% Reference: Akman, O.E., Watterson, S., Parton, A., Binns, N., Millar, A.J. & Ghazal, P. J. Roy. Soc. Interface, 9(74) (2012).
%
% DEPENDENCIES 
%
% None.
% 
% -------------------------------------------------------------------------
%
% Written by Ozgur Akman, University of Exeter, 2022
% O.E.Akman@exeter.ac.uk
%
% Part of the BDEtools package, © Akman Laboratory of Automated Biotechnology, 2021
%

yNow = [applygate1(Z(2, 2), g(2)) & Z(4, 6);
        applygate2(applygate1(Z(1, 1), g(1)), applygate1(Z(3, 5), g(5)), g(7));       
        applygate2(applygate1(Z(1, 3), g(3)), applygate1(Z(2, 4), g(4)), g(6)) & (Z(5, 7) | Z(6, 8));
        Z(4, 9);
        Z(5, 9);
        Z(6, 9)];

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