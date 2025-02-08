clear; close all; clc;

%% Define SNR range and Monte-Carlo parameters
snr_dB_range = 5:5:50;  % SNR in dB
num_MC = 5;            % number of Monte-Carlo trials per SNR

% Preallocate arrays for BER results.
BER_1x4_mat = zeros(num_MC, length(snr_dB_range));
BER_2x2_mat = zeros(num_MC, length(snr_dB_range));

%% Loop over SNR values
for snr_idx = 1:length(snr_dB_range)
    snr_val = snr_dB_range(snr_idx);
    fprintf('SNR = %d dB\n', snr_val);
    for mc = 1:num_MC
        %% --- Transmitter ---
        % Run your MIMO_OFDM_TX.m script (or function) to generate the TX signals.
        % It should define (at minimum):
        %   - tx_data, the transmitted bit stream
        %   - tx_vec_air_A and tx_vec_air_B for the 2×2 case.
        %   - For the 1×4 case, only tx_vec_air_A is used (the same waveform is sent on all 4 antennas,
        %     but possibly with different power scaling).
        %
        % Here we assume that MIMO_OFDM_TX.m is written as a script that uses a variable 'snr'
        % (in dB) to compute the noise power. We override it here.
        snr = snr_val;  % override TX snr
        MIMO_OFDM_TX;  % This script must generate the following variables:
        %   tx_data, tx_vec_air_A, tx_vec_air_B,
        %   and (for example) the downsampling filter "interp_filt2" and other parameters.
        
        %% --- Channel Simulation ---
        % (Note: The provided TX code simulates two cases.)
        %
        % For 1×4, we assume four received waveforms (rx_vec_air_1A, 1B, 1C, 1D) are generated as:
        rx_vec_air_1A = tx_vec_air_A;
        rx_vec_air_1B = 0.5 * tx_vec_air_A;
        rx_vec_air_1C = tx_vec_air_A;
        rx_vec_air_1D = 0.5 * tx_vec_air_A;
        % Append some zeros to mimic extra delay:
        rx_vec_air_1A = [rx_vec_air_1A, zeros(1,100)];
        rx_vec_air_1B = [rx_vec_air_1B, zeros(1,100)];
        rx_vec_air_1C = [rx_vec_air_1C, zeros(1,100)];
        rx_vec_air_1D = [rx_vec_air_1D, zeros(1,100)];
        
        % For 2×2, the TX signals are combined as:
        rx_vec_air_2A = tx_vec_air_A + 0.5 * tx_vec_air_B;
        rx_vec_air_2B = tx_vec_air_B + 0.5 * tx_vec_air_A;
        rx_vec_air_2A = [rx_vec_air_2A, zeros(1,100)];
        rx_vec_air_2B = [rx_vec_air_2B, zeros(1,100)];
        
        % Compute noise power based on the variance of the transmitted signal and desired SNR.
        % (Adjust the scaling as needed.)
        noise_power_1 = var(rx_vec_air_1A) * 10^(-snr_val/20);
        noise_power_2 = var(rx_vec_air_2A) * 10^(-snr_val/20);
        
        % Add AWGN to each received waveform.
        rx_vec_air_1A = rx_vec_air_1A + noise_power_1* (randn(size(rx_vec_air_1A)) + 1i*randn(size(rx_vec_air_1A)));
        rx_vec_air_1B = rx_vec_air_1B + noise_power_1* (randn(size(rx_vec_air_1B)) + 1i*randn(size(rx_vec_air_1B)));
        rx_vec_air_1C = rx_vec_air_1C + noise_power_1* (randn(size(rx_vec_air_1C)) + 1i*randn(size(rx_vec_air_1C)));
        rx_vec_air_1D = rx_vec_air_1D + noise_power_1* (randn(size(rx_vec_air_1D)) + 1i*randn(size(rx_vec_air_1D)));
        
        rx_vec_air_2A = rx_vec_air_2A + noise_power_2* (randn(size(rx_vec_air_2A)) + 1i*randn(size(rx_vec_air_2A)));
        rx_vec_air_2B = rx_vec_air_2B + noise_power_2* (randn(size(rx_vec_air_2B)) + 1i*randn(size(rx_vec_air_2B)));
        
        %% --- Downsampling ---
        % Use your generate_downsample_rx() function to downsample the received waveforms.
        rx_down_1A = generate_downsample_rx(rx_vec_air_1A, INTERP_RATE, interp_filt2);
        rx_down_1B = generate_downsample_rx(rx_vec_air_1B, INTERP_RATE, interp_filt2);
        rx_down_1C = generate_downsample_rx(rx_vec_air_1C, INTERP_RATE, interp_filt2);
        rx_down_1D = generate_downsample_rx(rx_vec_air_1D, INTERP_RATE, interp_filt2);
        
        rx_down_2A = generate_downsample_rx(rx_vec_air_2A, INTERP_RATE, interp_filt2);
        rx_down_2B = generate_downsample_rx(rx_vec_air_2B, INTERP_RATE, interp_filt2);
        
        %% --- Receiver ---
        % Call the completed MIMO receiver which processes both cases.
        % [decoded_data_1x4, decoded_data_2x2, BER_1x4, BER_2x2] = MIMO_OFDM_RX_test(...
        %     rx_down_1A, rx_down_1B, rx_down_1C, rx_down_1D, ...
        %     rx_down_2A, rx_down_2B, tx_data);
        MIMO_OFDM_RX_test;
        % Save the BER values for this Monte Carlo trial.
        BER_1x4_mat(mc, snr_idx) = ber_1x4;
        BER_2x2_mat(mc, snr_idx) = ber_2x2;
    end
end

%% Compute average BER and variance (across Monte Carlo trials) for each SNR
avg_BER_1x4 = mean(BER_1x4_mat,1);
var_BER_1x4 = var(BER_1x4_mat,0,1);
avg_BER_2x2 = mean(BER_2x2_mat,1);
var_BER_2x2 = var(BER_2x2_mat,0,1);

%% --- Plotting ---
figure;
semilogy(snr_dB_range, avg_BER_1x4, '-o', 'LineWidth', 2);
hold on;
semilogy(snr_dB_range, avg_BER_2x2, '-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Average BER');
legend('1×4 MRC', '2×2 ZF');
title('Average BER vs. SNR');

figure;
plot(snr_dB_range, var_BER_1x4, '-o', 'LineWidth', 2);
hold on;
plot(snr_dB_range, var_BER_2x2, '-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER Variance');
legend('1×4 MRC', '2×2 ZF');
title('BER Variance vs. SNR');
