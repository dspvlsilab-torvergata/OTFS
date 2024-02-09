r = zeros(N * M, 1);

for q = 1 : N * M
    for l = 1 : N
        if q >= l
            r(q) = r(q) + gs(l, q) * s(q - l + 1);
        end
    end
end