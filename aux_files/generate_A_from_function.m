% Copyright (c) 2016, Ryan Collingham.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Samples a function f(w) at discrete steps of the sample rate between w min
% and w max, and uses this to generate an approximate mapping of delta-like
% peaks of width dw. From here, an A matrix corresponding to this power-
% spectrum is generated for use in GLE simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = generate_A_from_function(f, dw, w_min, w_max, sample_rate)

% First, sample function between limits at dw.
[freq_bins, w] = sample_function(f, sample_rate, w_min, w_max);

% Re-scale the frequency bins to give correct height deltas.
freq_bins = freq_bins ./ scale_height(dw, w);

% Generate an A matrix using these freq bins
A = generate_A_from_frequencies_multiple_gamma(w, dw, freq_bins, 1);
end

function [freq_bins, w] = sample_function(f, dw, w_min, w_max)

% Sample function f(w) at in range w min - w max at each dw.
w = w_min:dw:w_max;
freq_bins = f(w);
end

function factor = scale_height(dw, w_0)
% Calculate the height of a delta-like peak of width dw at central freq
% w 0. Use this to scale the height of the peaks.
factor = 2 * sqrt(2 / pi) / dw * (dw^2 + 2 * w_0 .^2) ./ ...
    (dw^2 + 4 * w_0 .^ 2);
end