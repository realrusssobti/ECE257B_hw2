%% MIMO_OFDM_RX_TEST.m

%% Correlate for LTS
LTS_CORR_THRESH=.8;
DO_APPLY_CFO_CORRECTION=0;
DO_APPLY_SFO_CORRECTION=1;
DO_APPLY_PHASE_ERR_CORRECTION=1;

% For simplicity, we'll only use RFA for LTS correlation and peak
% discovery. A straightforward addition would be to repeat this process for
% RFB and combine the results for detection diversity.

% Complex cross correlation of Rx waveform with time-domain LTS
lts_corr = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_A)));

% Skip early and late samples - avoids occasional false positives from pre-AGC samples
lts_corr = lts_corr(32:end-32);

% Find all correlation peaks
lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));

% Select best candidate correlation peak as LTS-payload boundary
% In this MIMO example, we actually have 3 LTS symbols sent in a row.
% The first two are sent by RFA on the TX node and the last one was sent
% by RFB. We will actually look for the separation between the first and the
% last for synchronizing our starting index.

[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
[lts_last_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

% Stop if no valid correlation peak was found
if(isempty(lts_last_peak_index))
    fprintf('No LTS Correlation Peaks Found!\n');
    return;
end

% Set the sample indices of the payload symbols and preamble
% The "+32" here corresponds to the 32-sample cyclic prefix on the preamble LTS
% The "+192" corresponds to the length of the extra training symbols for MIMO channel estimation
mimo_training_ind = lts_peaks(max(lts_last_peak_index)) + 32;
payload_ind = mimo_training_ind + 192;

% Subtract of 2 full LTS sequences and one cyclic prefixes
% The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)
lts_ind = mimo_training_ind-160;

if(DO_APPLY_CFO_CORRECTION)
    %Extract LTS (not yet CFO corrected)
    rx_lts = raw_rx_dec_A(lts_ind : lts_ind+159); %Extract the first two LTS for CFO
    rx_lts1 = rx_lts(-64 + [97:160]);
    rx_lts2 = rx_lts( [97:160]);

    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts2 .* conj(rx_lts1))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

% Apply CFO correction to raw Rx waveforms
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec_A)-1]);
rx_dec_cfo_corr_A = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_B = raw_rx_dec_B .* rx_cfo_corr_t;


% MIMO Channel Estimatation
lts_ind_TXA_start = mimo_training_ind + 32 ;
lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;

lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 ;
lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;

rx_lts_AA = rx_dec_cfo_corr_A( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BA = rx_dec_cfo_corr_A( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AB = rx_dec_cfo_corr_B( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BB = rx_dec_cfo_corr_B( lts_ind_TXB_start:lts_ind_TXB_end );

rx_lts_AA_f = fft(rx_lts_AA);
rx_lts_BA_f = fft(rx_lts_BA);

rx_lts_AB_f = fft(rx_lts_AB);
rx_lts_BB_f = fft(rx_lts_BB);

%% Perform Channel estimation 

%% 1x1 mode 
rx_dec_cfo_corr_1A = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_1B = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_1C = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_1D = raw_rx_dec_A .* rx_cfo_corr_t;
% --- Channel Estimation using the legacy LTS ---
rx_lts_1A = rx_dec_cfo_corr_1A(lts_ind : lts_ind+159);
rx_lts_1B = rx_dec_cfo_corr_1B(lts_ind : lts_ind+159);
rx_lts_1C = rx_dec_cfo_corr_1C(lts_ind : lts_ind+159);
rx_lts_1D = rx_dec_cfo_corr_1D(lts_ind : lts_ind+159);
% Use the second LTS copy (samples 65:128)
rx_lts_1A_est = rx_lts_1A(65:128);
rx_lts_1B_est = rx_lts_1B(65:128);
rx_lts_1C_est = rx_lts_1C(65:128);
rx_lts_1D_est = rx_lts_1D(65:128);
lts_freq = fft(lts_t);
H_est_1A = fft(rx_lts_1A_est) ./ lts_freq;
H_est_1B = fft(rx_lts_1B_est) ./ lts_freq;
H_est_1C = fft(rx_lts_1C_est) ./ lts_freq;
H_est_1D = fft(rx_lts_1D_est) ./ lts_freq;

% --- Payload Extraction for 1×4 ---
payload_length = (CP_LEN + N_SC)*N_OFDM_SYMS;

payload_A_1 = rx_dec_cfo_corr_1A(payload_ind : payload_ind+payload_length-1);
payload_B_1 = rx_dec_cfo_corr_1B(payload_ind : payload_ind+payload_length-1);
payload_C_1 = rx_dec_cfo_corr_1C(payload_ind : payload_ind+payload_length-1);
payload_D_1 = rx_dec_cfo_corr_1D(payload_ind : payload_ind+payload_length-1);
payload_mat_1A = reshape(payload_A_1, CP_LEN+N_SC, N_OFDM_SYMS);
payload_mat_1B = reshape(payload_B_1, CP_LEN+N_SC, N_OFDM_SYMS);
payload_mat_1C = reshape(payload_C_1, CP_LEN+N_SC, N_OFDM_SYMS);
payload_mat_1D = reshape(payload_D_1, CP_LEN+N_SC, N_OFDM_SYMS);
% Remove CP and take FFT:
syms_f_mat_1A = fft(payload_mat_1A(CP_LEN+1:end, :), N_SC, 1);
syms_f_mat_1B = fft(payload_mat_1B(CP_LEN+1:end, :), N_SC, 1);
syms_f_mat_1C = fft(payload_mat_1C(CP_LEN+1:end, :), N_SC, 1);
syms_f_mat_1D = fft(payload_mat_1D(CP_LEN+1:end, :), N_SC, 1);

% --- Maximal-Ratio Combining (MRC) ---
rx_syms_case_1 = zeros(length(SC_IND_DATA)*N_OFDM_SYMS, 1);
for n = 1:N_OFDM_SYMS
    for idx = 1:length(SC_IND_DATA)
        k = SC_IND_DATA(idx);
        comb_num = conj(H_est_1A(k))*syms_f_mat_1A(k,n) + ...
                   conj(H_est_1B(k))*syms_f_mat_1B(k,n) + ...
                   conj(H_est_1C(k))*syms_f_mat_1C(k,n) + ...
                   conj(H_est_1D(k))*syms_f_mat_1D(k,n);
        comb_den = abs(H_est_1A(k))^2 + abs(H_est_1B(k))^2 + abs(H_est_1C(k))^2 + abs(H_est_1D(k))^2;
        rx_syms_case_1((n-1)*length(SC_IND_DATA) + idx) = comb_num / comb_den;
    end
end

%
%% 2x2 mode
% (Extract two 64-sample segments for TX A and TX B from rx_dec_cfo_corr_A and rx_dec_cfo_corr_B.)
lts_ind_TXA_start = mimo_training_ind + 32;
lts_ind_TXA_end   = lts_ind_TXA_start + 64 - 1;
lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32;
lts_ind_TXB_end   = lts_ind_TXB_start + 64 - 1;

rx_lts_AA = rx_dec_cfo_corr_A(lts_ind_TXA_start : lts_ind_TXA_end);
rx_lts_BA = rx_dec_cfo_corr_A(lts_ind_TXB_start : lts_ind_TXB_end);
rx_lts_AB = rx_dec_cfo_corr_B(lts_ind_TXA_start : lts_ind_TXA_end);
rx_lts_BB = rx_dec_cfo_corr_B(lts_ind_TXB_start : lts_ind_TXB_end);

% Compute FFT of the known LTS (assumed to be lts_t)
lts_freq = fft(lts_t);
H_est_TXA_rxA = fft(rx_lts_AA) ./ lts_freq;
H_est_TXB_rxA = fft(rx_lts_BA) ./ lts_freq;
H_est_TXA_rxB = fft(rx_lts_AB) ./ lts_freq;
H_est_TXB_rxB = fft(rx_lts_BB) ./ lts_freq;

% --- Payload Extraction ---
payload_length = (CP_LEN + N_SC)*N_OFDM_SYMS;
payload_A = rx_dec_cfo_corr_A(payload_ind : payload_ind+payload_length-1);
payload_B = rx_dec_cfo_corr_B(payload_ind : payload_ind+payload_length-1);
payload_mat_A = reshape(payload_A, CP_LEN+N_SC, N_OFDM_SYMS);
payload_mat_B = reshape(payload_B, CP_LEN+N_SC, N_OFDM_SYMS);
% Remove CP:
payload_mat_noCP_A = payload_mat_A(CP_LEN+1:end, :);
payload_mat_noCP_B = payload_mat_B(CP_LEN+1:end, :);
% FFT per OFDM symbol:
syms_f_mat_A = fft(payload_mat_noCP_A, N_SC, 1);
syms_f_mat_B = fft(payload_mat_noCP_B, N_SC, 1);

% --- Per-Subcarrier MIMO Equalization (zero-forcing) ---
% For each subcarrier (and each OFDM symbol) form the 2×2 channel matrix:
rx_syms_case_2 = zeros(length(SC_IND_DATA)*N_OFDM_SYMS, 2);
for n = 1:N_OFDM_SYMS
    for idx = 1:length(SC_IND_DATA)
        k = SC_IND_DATA(idx);
        y = [syms_f_mat_A(k, n); syms_f_mat_B(k, n)];
        H = [H_est_TXA_rxA(k), H_est_TXB_rxA(k); H_est_TXA_rxB(k), H_est_TXB_rxB(k)];
        s_hat = H \ y;  % zero-forcing equalization
        rx_syms_case_2((n-1)*length(SC_IND_DATA) + idx, :) = s_hat.';
    end
end
rx_syms_case_2_dec = rx_syms_case_2(:, 1);
rx_syms_case_1 = rx_syms_case_1(:,1);

% Because we only used Tx RFA to send pilots, we can do SISO equalization
% here. This is zero-forcing (just divide by chan estimates)
% syms_eq_mat_pilots = syms_f_mat_A ./ repmat(rx_H_est_AA.', 1, N_OFDM_SYMS);
syms_eq_mat_pilots = syms_f_mat_A;  % For pilot extraction

if DO_APPLY_SFO_CORRECTION
    % --- Step 1: Extract received pilot tones and estimate phase error ---
    rx_pilots = syms_eq_mat_pilots(SC_IND_PILOTS, :);
    % Remove known modulation by multiplying with the conjugate of the transmitted pilots
    pilot_phase_est = angle(rx_pilots .* conj(repmat(pilots_A, 1, N_OFDM_SYMS)));
    % Preallocate an SFO correction phase matrix (one phase value per subcarrier per OFDM symbol)
    pilot_phase_sfo_corr = zeros(N_SC, N_OFDM_SYMS);
    % For each OFDM symbol, interpolate the phase error from the pilot subcarrier indices to all subcarriers.
    for sym_idx = 1:N_OFDM_SYMS
        pilot_phase = pilot_phase_est(:, sym_idx);
        x  = SC_IND_PILOTS;
        xi = 1:N_SC;
        pilot_phase_sfo_corr(:, sym_idx) = interp1(x, pilot_phase, xi, 'linear', 'extrap');
    end
    % Form the SFO correction factor
    sfo_exp_matrix = exp(-1j * pilot_phase_sfo_corr);
    % Apply SFO correction to the entire frequency-domain symbol matrix for Rx A
    syms_f_mat_A = syms_eq_mat_pilots .* sfo_exp_matrix;
else
    pilot_phase_sfo_corr = zeros(N_SC, N_OFDM_SYMS);
end
%*This is optional* 
% Extract the pilots and calculate per-symbol phase error
if DO_APPLY_PHASE_ERR_CORRECTION
    pilots_f_mat = syms_eq_mat_pilots(SC_IND_PILOTS, :);
    pilot_phase_err = angle(mean(pilots_f_mat.*pilots_A));
else
	% Define an empty phase correction vector (used by plotting code below)
    pilot_phase_err = zeros(1, N_OFDM_SYMS);
end
pilot_phase_corr = repmat(exp(-1i*pilot_phase_err), N_SC, 1);

% Apply pilot phase correction to both received streams
syms_f_mat_pc_A = syms_f_mat_A .* pilot_phase_corr;
syms_f_mat_pc_B = syms_f_mat_B .* pilot_phase_corr;

% Perform combining for MIMO 1X4 and 2X2 
% you need to apply the MIMO equalization to each subcarrier separately and then perform combining

payload_syms_mat_A = rx_syms_case_2_dec(SC_IND_DATA, :);
payload_syms_mat_B = rx_syms_case_1(SC_IND_DATA, :);

%% perform demodulate or demapping post combined symbols 

% Select odd and even indices
odd_inds  = 1:2:length(rx_syms_case_1);
even_inds = 2:2:length(rx_syms_case_1);

% For odd indices: flip the sign of the real part only.
rx_syms_case_1(odd_inds) = complex(-real(rx_syms_case_1(odd_inds)), imag(rx_syms_case_1(odd_inds)));

% For even indices: take the complex conjugate.
rx_syms_case_1(even_inds) = conj(rx_syms_case_1(even_inds));
% rx_syms_case_1 = conj(rx_syms_case_1);

maxVal_1 = max(abs(rx_syms_case_2_dec(:)));
maxVal_2 = max(abs(rx_syms_case_1(:)));
% normalize case 1 ( this is a workaround )
    rx_syms_case_1 = rx_syms_case_1 / (maxVal_2 / maxVal_1);
    % rx_syms_case_1 = -1.*rx_syms_case_1;



% plot the demodulated output rx_syms_case_1 and rx_syms_case_2
figure(4);
scatter(real(rx_syms_case_1), imag(rx_syms_case_1),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

figure(5);
scatter(real(rx_syms_case_2), imag(rx_syms_case_2),'filled');
title(' Signal Space of received bits ( both of them )');
xlabel('I'); ylabel('Q');

figure(6);
scatter(real(rx_syms_case_2_dec), imag(rx_syms_case_2_dec),'filled');
title(' Signal Space of received bits ( both of them )');
xlabel('I'); ylabel('Q');


% FEC decoder for the rx_syms_case_1 and rx_syms_case_2

Demap_out_case_1 = demapper(rx_syms_case_1,MOD_ORDER,1);
Demap_out_case_2 = demapper(rx_syms_case_2_dec,MOD_ORDER,1);

% viterbi decoder
rx_data_final_1 = vitdec(Demap_out_case_1,trel,7,'trunc','hard');
rx_data_final_2 = vitdec(Demap_out_case_2,trel,7,'trunc','hard');

% rx_data is the final output corresponding to tx_data, which can be used
% to calculate BER

[number,ber_1x4] = biterr(tx_data_a,rx_data_final_1);
fprintf("1x4 bit error: %d \n", ber_1x4);
[number,ber_2x2] = biterr(tx_data_b,rx_data_final_2);
fprintf("2x2 bit error: %d \n", ber_2x2);
