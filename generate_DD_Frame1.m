function [X, GI] = generate_DD_Frame1(N, M, mod_size)

% number of information symbols in one frame:
N_syms_per_frame = N*M;

% number of information bits in one frame:
N_bits_per_frame = N_syms_per_frame * log2(mod_size);

% generate random bits:
tx_info_bits = randi([0,1], N_bits_per_frame, 1);
% tx_info_bits = [0; 0; 0; 1; 1; 1; 1; 0; 0; 0; 0; 1; 1; 1; 1; 0];

% QAM modulation:
tx_info_symbols = qammod(tx_info_bits, mod_size, 'gray', 'InputType', 'bit');

% Generate the MxN OTFS delay-Doppler frame:
X = reshape(tx_info_symbols, N, M);

GI = zeros(N, M);
GI(N/2, M/2) = 10;

end