function [I, J, S] = non_null_index(H, N, M)
S1 = zeros(N * M, 1);

for w = 1 : N * M
    S1(w) = nnz(abs(H(w, :)));
end

S = max(S1);

H1 = H.';

k = find(abs(H));
k1 = find(abs(H1));

I = reshape(k, S, N * M).';
J = reshape(k1, S, N * M).'; % la differenza è che la applico a H trasposta

for a = 1 : N * M
    for b = 1 : S
        I(a, b) = I(a, b) - ((a - 1) * (N * M));
        J(a, b) = J(a, b) - ((a - 1) * (N * M));
    end
end

end