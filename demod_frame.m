r_dt = reshape(r, N, M);
r_dd = (1 / M) * fft(r_dt, M, 2);