clc
clear
close all

% number of Doppler bins (time slots):
M = 4;

% number of delay bins (subcarriers):
N = 2;

%modulation size:
mod_size = 4;

% FFT Matrix:
Fn = dftmtx(N);
Fm = dftmtx(M);

% Identity Matrix:
In = eye(N);
Im = eye(M);

% Generate the N x M OTFS Delay-Doppler frame:
[X, GI] = generate_DD_Frame1(N, M, mod_size);

% Channel Parameters:
taps = 2;
g = [1 0.4];
l = [0 1];
k = [0 -1];

% Delay - Doppler Channel Matrix (Equazione 2.12):
Hdd = zeros(N, M);
for t = 1 : length(g)
    Hdd(l(t) + 1, k(t) + (M / 2)) = g(t);
end

% Delay - Time Response Channel (Equazione 4.8):
gs = zeros(N, N * M);
Z = exp((2 * pi * 1j)/(N * M));
for a = 1 : N     
    for q = 1 : N * M
        for i = 1 : taps
            if a == l(i) + 1
                gs(a, q) = gs(a, q) + g(i) * Z ^ (k(i) * (q - 1 - l(i)));
            end
        end
    end
end

%% Matrix:
% Permutation Matrix:
P = zeros(N * M, N * M);
for a = 1 : M
    for b = 1 : N
        E = zeros(N, M);
        E(b, a) = 1;
        P((a - 1) * N + 1 : a * N, (b - 1) * M + 1 : b * M) = E;
    end
end

% Delay - Time Channel Matrix:
G = zeros(N * M, N * M);
for q = 1 : N * M
	for ell = 1 : N
		if(q >= ell)
			G(q, q - ell + 1) = gs(ell, q);
		end
	end
end

% Delay - Doppler Channel Matrix:
H = (1 / M) * kron(In, Fm) * (P' * G * P) * kron(In, Fm');

%% Time Analysis:
% Trasmission signal (Equazione 4.3):
GI_dt = M * ifft(GI, M, 2);        
s = reshape(GI_dt, N * M, 1);

X_dt = M * ifft(X, M , 2);
s1 = reshape(X_dt, N * M, 1);

%% Delay - Doppler Analysis:
% Delay - Doppler Response Data using the Delay - Doppler Channel Matrix H (Equazione 4.58):
GI_dd = (1 / M) * (kron(In, Fm) * (P' * s));
X_dd = (1 / M) * (kron(In, Fm) * (P' * s1));

Y_GI_dd = H * GI_dd;
Y_X_dd = H * X_dd; 

r_GI_dd = reshape(Y_GI_dd, M, N).';
r_X_dd = reshape(Y_X_dd, M, N).';

%% Message Passing:
[I, J, S] = non_null_index(H, N, M);

pmf = (1 / mod_size) * ones(N * M, N * M, mod_size); % sarebbe P(c,d)
pmf1 = ones(N * M, N * M, mod_size); % sarebbe P(c,d) ondina
pmf2 = ones(N * M, mod_size); % sarebbe P(c)

damping_factor = 0.5;

media = zeros(N * M, N * M);

varianza1 = 0;
varianza2 = 0;
varianza = zeros(N * M, N * M);

somma = 0;
somma1 = 0;

gamma = 0.5;
counter = 0;

epsilon = 0.2;
reg = 0;

bit = [0; 0; 1; 0; 1; 1; 0; 1];
simboli_mod = qammod(bit, mod_size, 'gray', 'InputType', 'bit');

index_possible = zeros(N * M, 1);
probability = zeros(N * M, 1);
X_tilde = zeros(N * M, 1);

% Numero di iterazioni:
for i = 1 : 10

for d = 1 : (N * M) 

for c = 1 : (N * M) 

% Calcolo della Media:
for e = 1 : S
    if I(d, e) ~= c
        for j = 1 : length(simboli_mod)
            media(d, c) = media(d, c) + pmf(c, d, j) * simboli_mod(j) * H(d, I(d, e));
        end
    end
end

% Calcolo della Varianza:
for e = 1 : S
    if I(d, e) ~= c
        for j = 1 : length(simboli_mod)
            varianza1 = varianza1 + pmf(c, d, j) * (abs(simboli_mod(j)) ^ 2) * (abs(H(d, I(d, e))) ^ 2);
            varianza2 = varianza2 + pmf(c, d, j) * simboli_mod(j) * H(d, I(d, e));
        end
    varianza2 = abs(varianza2) ^ 2;
    varianza(d, c) = varianza(d, c) + (varianza1 - varianza2);
    end
end

varianza1 = 0;
varianza2 = 0;

end % fine c

end % fine d

for t = 1 : N * M % c

for g = 1 : N * M % d

for f = 1 : length(simboli_mod)
    for m = 1 : S    
        if J(t, m) ~= g
           for n = 1 : length(simboli_mod) 
                somma = somma + exp(-(abs(Y_X_dd(J(t, m)) - media(J(t, m), t) - (H(J(t, m), t) * simboli_mod(n))) ^ 2) / varianza(J(t, m), t));     
           end
        pmf1(t, g, f) = pmf1(t, g, f) * (exp(-(abs(Y_X_dd(J(t, m)) - media(J(t, m), t) - (H(J(t, m), t) * simboli_mod(f))) ^ 2) / varianza(J(t, m), t))) / (somma);
        somma = 0;
        end
    end
    pmf(t, g, f) = pmf1(t, g, f) * damping_factor + (1 - damping_factor) * pmf(t, g, f);
end

end % fine d 

for r = 1 : length(simboli_mod)
    for p = 1 : S    
        for q = 1 : length(simboli_mod) 
            somma1 = somma1 + exp(-(abs(Y_X_dd(J(t, p)) - media(J(t, p), t) - (H(J(t, p), t) * simboli_mod(q))) ^ 2) / varianza(J(t, p), t));     
        end
    pmf2(t, r) = pmf2(t, r) * (exp(-(abs(Y_X_dd(J(t, p)) - media(J(t, p), t) - (H(J(t, p), t) * simboli_mod(r))) ^ 2) / varianza(J(t, p), t))) / (somma1);
    somma1 = 0;
    end
end

if max(pmf2(t, :)) > (1 - gamma)
    counter = counter + 1;
end

end % fine c

eta = (1 / (N * M)) * counter;

if reg < eta
    reg = eta;
end

for u = 1 : N * M
    [probability(u), index_possible(u)] = max(pmf2(u, :));
end

if eta == 1 || eta < reg - epsilon
     break
end

pmf1 = ones(N * M, N * M, mod_size); % sarebbe P(c,d) ondina
pmf2 = ones(N * M, mod_size); % sarebbe P(c)

end % fine iterazioni

for v = 1 : N * M
X_tilde(v) = simboli_mod(index_possible(v));
end

%% Plot:
% figure
% bar3(abs(Hdd));
% grid on
% ylabel('Delay');
% xlabel('Doppler');
% title('Delay Doppler Channel Matrix:');
% 
% figure
% subplot(2, 2, 1)
% bar3(abs(GI));
% grid on
% ylabel('Delay');
% xlabel('Doppler');
% title('Delay Doppler Input Pilot:');
% 
% subplot(2, 2, 2)
% bar3(abs(X));
% grid on
% ylabel('Delay');
% xlabel('Doppler');
% title('Delay Doppler Input Data:');
% 
% subplot(2, 2, 3)
% bar3(abs(r_GI_dd));
% grid on
% ylabel('Delay');
% xlabel('Doppler');
% title('Delay Doppler Output Pilot:');
% 
% subplot(2, 2, 4)
% bar3(abs(r_X_dd));
% grid on
% ylabel('Delay');
% xlabel('Doppler');
% title('Delay Doppler Output Data:');
