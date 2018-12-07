% Copyright (c) 2018, Nadav Avidor.
% Copyright (c) 2016, Ryan Collingham.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate an A matrix for a comb of delta-like peaks at frequencies w list,
% each with width dw and height proportional to each of gamma list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = generate_A_from_frequencies_multiple_gamma(w_list, dw_list, gamma_list, tau_list)

if length(w_list) > 1
    if length(dw_list) == 1, dw_list = repmat(dw_list,size(w_list)); end
    if length(gamma_list) == 1, gamma_list = repmat(gamma_list,size(w_list)); end
    if length(tau_list) == 1, tau_list = repmat(tau_list,size(w_list)); end
end

% Determine the size of A and allocate a matrix of zeros.
n = 2 * length(w_list) + 1;
A = zeros(n, n);

% Now fill out the rest of the matrix.
for i = 2:n
A = fill_in_values(A, i, gamma_list(floor(i / 2)), tau_list(floor(i / 2)), dw_list(floor(i / 2)), w_list(floor(i / 2)));
end

end

function A = fill_in_values(A, i, gamma, tau, dw, w)
% Fills out the i'th row, collumn and diagonal elements of A.
% First row
A(1, i) = sqrt(gamma / tau);
% First collumn
A(i, 1) = -sqrt(gamma / tau);
% Populate the rest of the main diagonal with dw
A(i, i) = dw;
if mod(i, 2) == 0
% Populate the next diagonal along with the frequencies contained
% in w list.
A(i + 1, i) = w;
% Finally, populate the preceding diagonal with the negative of the
% frequencies contained in w list.
A(i, i + 1) = - w;
end
end