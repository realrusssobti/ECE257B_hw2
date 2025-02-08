%% MyOfdmReceiver.m
function [decoded_data]= MyOfdmReceiver(rcvd_data) 
 %% run transmitter code to load sts and lts and other parameters 
 OFDM_TX; 
%% Rx processing params
%rx_data = rcvd_data;   % When passing the right data and used as a function
rx_data = raw_rx_data;          % run OFDM tx code to get raw_rx_dec
rx_data = rx_data ./ max(abs(rx_data));
LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation

% Usage: Find all peaks whose magnitude is greater than 0.8 times
USE_STS_PACKET_DETECTION=1;
SFO_ENABLE=1;
PIECEWISE_INTERP = 1;
USE_ONLYLTS_CFO=0;
USE_ONELTS_CHAN_EST=0;
USE_LTS_FFT=0;

% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet

%% Packet DetectionUSE_STS_PACKET_DETECTION=1;


% ideas: Cross correlation of received signal with LTS or use STS to detect the packet?

length_samples= length(rx_data) - 200;

% Store 0 in first 16 elements since those cant be autocorrelated
self_corr_sample=length(sts_t);
cross_corr_sample=0;
self_corr_output=zeros(size(rx_data));
sts_acr_output=zeros(size(rx_data));
lts_acr_output=zeros(size(rx_data));

while(self_corr_sample < length_samples)
    
    self_corr_output(self_corr_sample+1)= rx_data(self_corr_sample-length(sts_t) + (1:length(sts_t))) * rx_data(self_corr_sample + (1:length(sts_t)))' ./norm(rx_data(self_corr_sample+(1:length(sts_t))));
    sts_acr_output(cross_corr_sample+1)= sts_t * rx_data(cross_corr_sample + (1:length(sts_t)))' ./norm(rx_data(cross_corr_sample+(1:length(sts_t))));  
    lts_acr_output(cross_corr_sample+1)= lts_t * rx_data(cross_corr_sample + (1:length(lts_t)))' ./norm(rx_data(cross_corr_sample+(1:length(lts_t)))); 
    self_corr_sample= self_corr_sample+1;
    cross_corr_sample=cross_corr_sample+1;
end

self_corr_output= self_corr_output./max(abs(self_corr_output));
sts_acr_output= sts_acr_output./max(abs(sts_acr_output));
lts_acr_output= lts_acr_output./max(abs(lts_acr_output));

figure;
plot((abs(self_corr_output))); 
hold on 
plot(abs(rx_data)); 
hold off

% 30 Peaks observed
figure;
sts_peaks =  find(abs(sts_acr_output)>LTS_CORR_THRESH);
plot(abs(sts_acr_output));

num_sts = 30;
sts_first_peak = 1;
if(numel(sts_peaks)>=num_sts)
    sts_sanity =0;
    while (sts_first_peak<=numel(sts_peaks))
        if((sts_first_peak+num_sts-1)>numel(sts_peaks)) 
            break
        end
        sts_candidates = sts_peaks(sts_first_peak+(0:num_sts-1));
        sts_peak_diffs = sts_candidates - circshift(sts_candidates,1);
        sts_peak_diffs = abs(sts_peak_diffs(2:end)-length(sts_t));
        decision_vec=size(find(sts_peak_diffs>0));
        if(decision_vec(2)==0)
            sts_sanity = 1;
            break
        else
            sts_first_peak=sts_first_peak+1;
        end
    end
    if(sts_sanity==0)
        sprintf("Packet Detection unable to find 30 STS peaks spaced 16")
        return
    end
else
    numel(sts_peaks)
    sprintf("Packet Detection unable to find 30 STS peaks")
    return
end

% 2 Peaks observed
figure;
lts_peaks =  find(abs(lts_acr_output)>LTS_CORR_THRESH);
plot(abs(lts_acr_output));
num_lts = 2;
lts_first_peak = 1;
if(numel(lts_peaks)>=num_lts)
    lts_sanity =0;
    while (lts_first_peak<=numel(lts_peaks))
        if((lts_first_peak+num_lts-1)>numel(lts_peaks)) 
            break
        end
        lts_candidates = lts_peaks(lts_first_peak+(0:num_lts-1));
        lts_peak_diffs = lts_candidates - circshift(lts_candidates,1);
        lts_peak_diffs = abs(lts_peak_diffs(2:end)-length(lts_t));
        decision_vec=size(find(lts_peak_diffs>0));
        if(decision_vec(2)==0)
            lts_sanity = 1;
            break
        else
            lts_first_peak=lts_first_peak+1;
        end
    end
    if(lts_sanity==0)
        sprintf("Packet Detection unable to find 2 LTS peaks spaced 64")
        return
    end
else
    sprintf("Packet Detection unable to find 2 LTS peaks spaced 64")
    return
end


if(USE_STS_PACKET_DETECTION==1)
    packet_start =  30*length(sts_t)+sts_peaks(sts_first_peak)+2.5*length(lts_t);
    lts_start_index = 30*length(sts_t)+sts_peaks(sts_first_peak);
else
    packet_start =  2*length(lts_t)+lts_peaks(lts_first_peak);
    lts_start_index = lts_peaks(lts_first_peak)-0.5*length(lts_t);
end

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices
% Find peaks

%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)

if(USE_ONLYLTS_CFO)
    lts_len = length(lts_t);
    rx_lts = rx_data(lts_start_index:lts_start_index+2.5*lts_len);
    % Remove the end parts of lts to avoid occasional false positives from pre-AGC samples
    rx_lts_end_removed = rx_lts(0.25*lts_len:2.25*lts_len);
    cfo_est = mean((unwrap(angle(rx_lts_end_removed(1:lts_len)./...
        rx_lts_end_removed(lts_len+1:2*lts_len))))/(2*pi*lts_len));
    cfo_removed_sig = exp(1i*2*pi*(cfo_est)*(1:numel(rx_data))).*rx_data;
else
    sts_len = length(sts_t);
    rx_sts = rx_data(sts_peaks(sts_first_peak):sts_peaks(sts_first_peak)+30*sts_len-1);
    % Remove the end parts of sts to avoid occasional false positives from pre-AGC samples
    rx_sts_end_removed = rx_sts(sts_len:29*sts_len-1);
    rx_sts_corr_arr = rx_sts(2*sts_len:30*sts_len-1);
    coarse_cfo_est = mean((unwrap(angle(rx_sts_corr_arr'.*rx_sts_end_removed.')))/(2*pi*sts_len));
    coarse_cfo_removed_sig = exp(1i*2*pi*(coarse_cfo_est)*(1:numel(rx_data))).*rx_data;
    
    lts_len = length(lts_t);
    rx_lts = coarse_cfo_removed_sig(lts_start_index:lts_start_index+2.5*lts_len-1);
    % Remove the end parts of lts to avoid occasional false positives from pre-AGC samples
    rx_lts_end_removed = rx_lts(0.25*lts_len:2.25*lts_len);
    rx_lts_1 = rx_lts_end_removed(1:lts_len);
    rx_lts_2 = rx_lts_end_removed(lts_len+1:2*lts_len);
    cfo_est = mean((unwrap(angle(rx_lts_2'.*rx_lts_1.')))/(2*pi*lts_len));
    cfo_removed_sig = exp(1i*2*pi*(cfo_est)*(1:numel(rx_data))).*coarse_cfo_removed_sig;
end
% Output: Packet with each value multiplied by CFO correction factor

%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help
rx_payload = cfo_removed_sig(packet_start:packet_start+(CP_LEN+N_SC)*N_OFDM_SYMS-1);
rx_payload_mat = reshape(rx_payload,[CP_LEN+N_SC,N_OFDM_SYMS]);
rx_payload_mat_cp_removed = rx_payload_mat(17:end,:);

% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)


%% FFT
% Refer to IFFT perfomed at TX
rx_symbols = fft(rx_payload_mat_cp_removed,N_SC,1);

% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix 
% Output : Symbol equalized matrix in frequency domain of same size
lts_cfo_cp_removed = cfo_removed_sig(lts_start_index+0.5*lts_len:packet_start-1);
lts_mat = reshape(lts_cfo_cp_removed,[lts_len,2]);
rx_lts_f = fft(lts_mat,N_SC,1);
non_zero_lts_freqs = find(lts_f~=0);
chan_lts_est = ones(N_SC,1);

if(USE_ONELTS_CHAN_EST)
    if(USE_LTS_FFT)
        chan_lts_est = mean(rx_lts_f(:,1)./fft(lts_t)',2);
    else
        chan_lts_est(non_zero_lts_freqs) = mean(rx_lts_f(non_zero_lts_freqs,1)./lts_f(non_zero_lts_freqs)',2);
    end
else
    if(USE_LTS_FFT)
        chan_lts_est = mean(rx_lts_f./fft(lts_t)',2);
    else
        chan_lts_est(non_zero_lts_freqs) = mean(rx_lts_f(non_zero_lts_freqs,:)./lts_f(non_zero_lts_freqs)',2);
    end
end
equalized_symbs = rx_symbols./chan_lts_est;
%% Advanced topics: 
%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.
% Output: Symbol equalized matrix with pilot phase correction applied

equalized_pilots = equalized_symbs(SC_IND_PILOTS,:);
pilot_phases = angle(equalized_pilots./pilots);
sfo_phase_matrix = ones(N_SC,N_OFDM_SYMS);

if(PIECEWISE_INTERP)
    %Interpolation 1
    interp_start = [0 pilot_phases(4,1:end-1)];
    interp_end = pilot_phases(1,:);
    x = [SC_IND_PILOTS(4)-N_SC-1,SC_IND_PILOTS(1)]';
    v = [interp_start; interp_end]; 
    xq = 1:SC_IND_PILOTS(1);
    sfo_phase_matrix(xq,:) = interp1(x,v,xq);
    %Interpolation 2
    interp_start = pilot_phases(1,:);
    interp_end = pilot_phases(2,:);
    x = [SC_IND_PILOTS(1),SC_IND_PILOTS(2)]';
    v = [interp_start; interp_end]; 
    xq = SC_IND_PILOTS(1)+1:SC_IND_PILOTS(2);
    sfo_phase_matrix(xq,:) = interp1(x,v,xq);
    %Interpolation 3
    interp_start = pilot_phases(2,:);
    interp_end = pilot_phases(3,:);
    x = [SC_IND_PILOTS(2),SC_IND_PILOTS(3)]';
    v = [interp_start; interp_end]; 
    xq = SC_IND_PILOTS(2)+1:SC_IND_PILOTS(3);
    sfo_phase_matrix(xq,:) = interp1(x,v,xq);
    %Interpolation 4
    interp_start = pilot_phases(3,:);
    interp_end = pilot_phases(4,:);
    x = [SC_IND_PILOTS(3),SC_IND_PILOTS(4)]';
    v = [interp_start; interp_end]; 
    xq = SC_IND_PILOTS(3)+1:SC_IND_PILOTS(4);
    sfo_phase_matrix(xq,:) = interp1(x,v,xq);
    %Interpolation 5
    interp_start = pilot_phases(4,:);
    interp_end = [pilot_phases(1,2:end) 0];
    x = [SC_IND_PILOTS(4),SC_IND_PILOTS(1)+N_SC]';
    v = [interp_start; interp_end]; 
    xq = SC_IND_PILOTS(4)+1:N_SC;
    sfo_phase_matrix(xq,:) = interp1(x,v,xq);
else
    for i=(1:N_OFDM_SYMS)
        interp_x = [8,22,-21,-7];
        interp_y = pilot_phases(:,i);
        P = polyfit(interp_x,interp_y',1);
        subC = (-32:31);
        sfo_phase_matrix(:,i) = P(1)*subC+P(2);
    end
end

sfo_exp_matrix = exp(-1j*sfo_phase_matrix);
sfo_corrected_sig = sfo_exp_matrix.*equalized_symbs;

%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error

% Output: Symbol equalized matrix with pilot phase correction applied

% Remove pilots and flatten the matrix to a vector rx_syms
%% Demodulation
if(SFO_ENABLE)
    corrected_datasig = sfo_corrected_sig(SC_IND_DATA,:);
else
    corrected_datasig = equalized_symbs(SC_IND_DATA,:);
end
rx_syms = corrected_datasig(:);
figure;
scatter(real(rx_syms), imag(rx_syms),'filled');
title(' Signal Space of received bits');
xlabel('I'); ylabel('Q');

% FEC decoder
Demap_out = demapper(rx_syms,MOD_ORDER,1);

% viterbi decoder
decoded_data = vitdec(Demap_out,trel,7,'trunc','hard');


% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER
BER = sum(abs(decoded_data-tx_data))/length(tx_data);
