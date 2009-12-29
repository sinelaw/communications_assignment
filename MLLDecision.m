function [ r_i, A_r ] = MLLDecision( q_k, Symbols )
%MLLDecision Takes a vector of pairs (2 x N) of values q_k = [[q1_ki; q1_kq], [q2_ki; q2_kq], ...] from receive filter and finds the nearest symbol.
%   A_r = symbol found, r_i = index of symbol found in Symbols vector

q = q_k(1,:) + 1j * q_k(2,:);
dist = zeros(length(Symbols),length(q));
for i = 1:length(Symbols)
    dist(i,:) = abs(q - Symbols(i)); % no need to calculate abs^2, minimum is the same for abs.
end

[min_val, min_index] = min(dist);

A_r = Symbols(min_index);
r_i = min_index;

end

