%Modulazione

X_dt = M * ifft(X, M, 2); %antitrasformata
s1 = reshape(X_dt, N * M, 1); %parallelo seriale