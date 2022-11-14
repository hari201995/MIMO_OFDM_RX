function [BER1x4,BER2x2] = MIMO_OFDM_RX_test(snr)

%% Correlate for LTS
[rx_vec_dec_1A,rx_vec_dec_1B,rx_vec_dec_1C,rx_vec_dec_1D,raw_rx_dec_2A,raw_rx_dec_2B,tx_data_a,pilots_A,MOD_ORDER] = MIMO_OFDM_TX(snr);
BER1x4=0;
BER2x2=0;
% Waveform params
N_OFDM_SYMS             = 500;         % Number of OFDM symbols
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

LTS_CORR_THRESH=.8;
DO_APPLY_CFO_CORRECTION=1;
DO_APPLY_SFO_CORRECTION=1;
DO_APPLY_PHASE_ERR_CORRECTION=1;
trel = poly2trellis(7, [171 133]);              % Define trellis
mimo_1x4 =1;
mimo_2x2 =1;

lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

% For simplicity, we'll only use RFA for LTS correlation and peak
% discovery. A straightforward addition would be to repeat this process for
% RFB and combine the results for detection diversity.

% Complex cross correlation of Rx waveform with time-domain LTS

                                                    %1x4 MIMO

lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx_vec_dec_1A)));

                                                    %2x2 MIMO

lts_corr_rx1 = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_2A)));
lts_corr_rx2 = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_2B)));

% Skip early and late samples - avoids occasional false positives from pre-AGC samples

lts_corr = lts_corr(32:end-32);                                                                           % 1x4 MIMO %
lts_corr_2x2 = lts_corr_rx1(32:end-32);                                                            % 2x2 MIMO %


% Find all correlation peaks

lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));                  % 1x4 MIMO %
lts_peaks_2x2 = find(lts_corr > LTS_CORR_THRESH*max(lts_corr_2x2));  % 2x2 MIMO %
 
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
    
                                                            %1x4 MIMO%
mimo_training_ind = lts_peaks(max(lts_last_peak_index)) + 31;
payload_ind = mimo_training_ind + 192;

                                                               %2x2 MIMO%
mimo_training_ind_2 = lts_peaks(max(lts_last_peak_index)) + 31;
payload_ind_2= mimo_training_ind + 192;                

% Subtract of 2 full LTS sequences and one cyclic prefixes
% The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)

lts_ind = mimo_training_ind-160;                                                                    % 1x4 MIMO %
lts_ind_2x2 = mimo_training_ind-160;                                                            % 2x2 MIMO %

                                                            %% CFO Correction %%
                                                                       % 1x4 MIMO % 

if (DO_APPLY_CFO_CORRECTION) && (mimo_1x4)
    %Extract LTS (not yet CFO corrected)
    rx_lts = rx_vec_dec_1A(lts_ind : lts_ind+159); %Extract the first two LTS for CFO
    rx_lts1 = rx_lts(-64 + [97:160]);
    rx_lts2 = rx_lts( [97:160]);

    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts2 .* conj(rx_lts1))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

if (DO_APPLY_CFO_CORRECTION) && (mimo_2x2)
    %Extract LTS (not yet CFO corrected)
    rx_lts = rx_vec_dec_1A(lts_ind : lts_ind+159); %Extract the first two LTS for CFO
    rx_lts1 = rx_lts(-64 + [97:160]);
    rx_lts2 = rx_lts( [97:160]);

    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts2 .* conj(rx_lts1))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

% Apply CFO correction to raw Rx waveforms
                                         
                                            %1x4 Mimo
if mimo_1x4
     rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(rx_vec_dec_1A)-1]);
    rx_dec_cfo_corr_1A = rx_vec_dec_1A .* rx_cfo_corr_t;
    rx_dec_cfo_corr_1B = rx_vec_dec_1B.* rx_cfo_corr_t;
    rx_dec_cfo_corr_1C = rx_vec_dec_1C .* rx_cfo_corr_t;
    rx_dec_cfo_corr_1D = rx_vec_dec_1D .* rx_cfo_corr_t;
end
                                            %2x2 Mimo
 if mimo_2x2
         rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(rx_vec_dec_1A)-1]);
         rx_dec_cfo_corr_2A =  raw_rx_dec_2A .* rx_cfo_corr_t;
         rx_dec_cfo_corr_2B =  raw_rx_dec_2B .* rx_cfo_corr_t;
 end
                              %% MIMO Channel Estimatation %%
                                            %% 1x4 MIMO %%
 if mimo_1x4
         lts_ind_1_start= mimo_training_ind - 128+1;
         lts_ind_1_end = lts_ind_1_start +64-1;
        
         lts_ind_2_start = mimo_training_ind - 64+1;
         lts_ind_2_end = lts_ind_2_start +64-1;
         
         lts_ind_3_start= mimo_training_ind +32+1;
         lts_ind_3_end = lts_ind_3_start +64-1;
        
                                           % 2x2 MIMO %
        lts_ind_TXA_start = mimo_training_ind + 32 ;
        lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;
        
        lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 ;
        lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;
                                       
                                            % 1x4 Channel estimation 
        % LTF extraction
        rx_lts_1_A1A =  rx_dec_cfo_corr_1A( lts_ind_1_start:lts_ind_1_end );
        rx_lts_2_A1A =  rx_dec_cfo_corr_1A( lts_ind_2_start:lts_ind_2_end );
        rx_lts_3_A1A =  rx_dec_cfo_corr_1A( lts_ind_3_start:lts_ind_3_end );
        
        rx_lts_1_A1B =  rx_dec_cfo_corr_1B( lts_ind_1_start:lts_ind_1_end );
        rx_lts_2_A1B =  rx_dec_cfo_corr_1B( lts_ind_2_start:lts_ind_2_end );
        rx_lts_3_A1B =  rx_dec_cfo_corr_1B( lts_ind_3_start:lts_ind_3_end );
        
        rx_lts_1_A1C =  rx_dec_cfo_corr_1C( lts_ind_1_start:lts_ind_1_end );
        rx_lts_2_A1C =  rx_dec_cfo_corr_1C( lts_ind_2_start:lts_ind_2_end );
        rx_lts_3_A1C =  rx_dec_cfo_corr_1C( lts_ind_3_start:lts_ind_3_end );
        
        rx_lts_1_A1D =  rx_dec_cfo_corr_1D( lts_ind_1_start:lts_ind_1_end );
        rx_lts_2_A1D =  rx_dec_cfo_corr_1D( lts_ind_2_start:lts_ind_2_end );
        rx_lts_3_A1D =  rx_dec_cfo_corr_1D( lts_ind_3_start:lts_ind_3_end );
        
        rx_lts_1_A1A_f = fft(rx_lts_1_A1A);
        rx_lts_2_A1A_f = fft(rx_lts_2_A1A);
        rx_lts_3_A1A_f = fft(rx_lts_3_A1A);
        
        rx_lts_1_A1B_f = fft(rx_lts_1_A1B);
        rx_lts_2_A1B_f = fft(rx_lts_2_A1B);
        rx_lts_3_A1B_f = fft(rx_lts_3_A1B);
        
        rx_lts_1_A1C_f = fft(rx_lts_1_A1C);
        rx_lts_2_A1C_f = fft(rx_lts_2_A1C);
        rx_lts_3_A1C_f = fft(rx_lts_3_A1C);
        
        rx_lts_1_A1D_f = fft(rx_lts_1_A1D);
        rx_lts_2_A1D_f = fft(rx_lts_2_A1D);
        rx_lts_3_A1D_f = fft(rx_lts_3_A1D);
 end

                                    %2x2 MIMO Channel Estimation.

if mimo_2x2 
        lts_1_2A_start = mimo_training_ind_2-128+1;
        lts_1_2A_end = lts_1_2A_start + 64 -1;
        
        lts_2_2A_start = mimo_training_ind_2-64+1;
        lts_2_2A_end = lts_2_2A_start + 64 -1;
        
        lts_3_2A_start = mimo_training_ind_2+32+1;
        lts_3_2A_end = lts_3_2A_start + 64 -1;
        
        lts_4_2A_start = mimo_training_ind_2+96+32+1;
        lts_4_2A_end = lts_4_2A_start + 64 -1;

        rx_samples_lts_1_2A = rx_dec_cfo_corr_2A( lts_1_2A_start:lts_1_2A_end );
        rx_samples_lts_2_2A = rx_dec_cfo_corr_2A( lts_2_2A_start:lts_2_2A_end );
        rx_samples_lts_3_2A = rx_dec_cfo_corr_2A( lts_3_2A_start:lts_3_2A_end );
        rx_samples_lts_4_2A = rx_dec_cfo_corr_2A( lts_4_2A_start:lts_4_2A_end );

        rx_samples_lts_1_2B = rx_dec_cfo_corr_2B( lts_1_2A_start:lts_1_2A_end );
        rx_samples_lts_2_2B = rx_dec_cfo_corr_2B( lts_2_2A_start:lts_2_2A_end );
        rx_samples_lts_3_2B = rx_dec_cfo_corr_2B( lts_3_2A_start:lts_3_2A_end );
        rx_samples_lts_4_2B = rx_dec_cfo_corr_2B( lts_4_2A_start:lts_4_2A_end );

        rx_lts_1_2A_f_domain = fft(rx_samples_lts_1_2A);
        rx_lts_2_2A_f_domain = fft( rx_samples_lts_2_2A);
        rx_lts_3_2A_f_domain = fft(rx_samples_lts_3_2A);
        rx_lts_4_2A_f_domain = fft(rx_samples_lts_4_2A);

         rx_lts_1_2B_f_domain = fft(rx_samples_lts_1_2B);
         rx_lts_2_2B_f_domain = fft(rx_samples_lts_2_2B);
         rx_lts_3_2B_f_domain = fft(rx_samples_lts_3_2B);
         rx_lts_4_2B_f_domain = fft(rx_samples_lts_4_2B);
        

end
                                                                   %% Perform Channel estimation  %%

                                                               % 1x4 MIMO % 
 if mimo_1x4

    H_cap_A1A_f  = 0.5*(rx_lts_2_A1A_f+rx_lts_3_A1A_f) ./ lts_f;
    H_cap_A1B_f  = 0.5*(rx_lts_2_A1B_f+rx_lts_3_A1B_f) ./ lts_f;
    H_cap_A1C_f  = 0.5*(rx_lts_2_A1C_f+rx_lts_3_A1C_f) ./ lts_f;
    H_cap_A1D_f  = 0.5*(rx_lts_2_A1D_f+rx_lts_3_A1D_f) ./ lts_f;
end
                                                                % 2x2 MIMO % 
if mimo_2x2

    H_cap_11_A_f = 0.5 * (rx_lts_1_2A_f_domain + rx_lts_2_2A_f_domain ) ./ lts_f;
    H_cap_12_A_f = rx_lts_4_2A_f_domain ./ lts_f;

    H_cap_22_B_f = rx_lts_4_2B_f_domain ./ lts_f;
    H_cap_21_B_f = 0.5*(rx_lts_1_2B_f_domain + rx_lts_2_2B_f_domain) ./ lts_f;

end
                                                                    %% Rx payload processing, Perform combining for 1X4 and 2X2 separately  

% Extract the payload samples (integral number of OFDM symbols following preamble)

if mimo_1x4
                                                                                        % 1x4 MIMO %

rx_payload_mat_with_cp_1A = reshape(rx_dec_cfo_corr_1A(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);
rx_payload_mat_with_cp_1B = reshape(rx_dec_cfo_corr_1B(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);
rx_payload_mat_with_cp_1C = reshape(rx_dec_cfo_corr_1C(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);
rx_payload_mat_with_cp_1D = reshape(rx_dec_cfo_corr_1D(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);

                                                                                          % Removing CP %
rx_payload_mat_without_cp_1A = rx_payload_mat_with_cp_1A(CP_LEN+[1:N_SC]:end,:);
rx_payload_mat_without_cp_1B = rx_payload_mat_with_cp_1A(CP_LEN+[1:N_SC]:end,:);
rx_payload_mat_without_cp_1C = rx_payload_mat_with_cp_1A(CP_LEN+[1:N_SC]:end,:);
rx_payload_mat_without_cp_1D = rx_payload_mat_with_cp_1A(CP_LEN+[1:N_SC]:end,:);
                                                                                           
                                                                                        % FFT 1x4 MIMO %

syms_f_mat_1A = fft(rx_payload_mat_without_cp_1A, N_SC, 1);
syms_f_mat_1B = fft(rx_payload_mat_without_cp_1B, N_SC, 1);
syms_f_mat_1C = fft(rx_payload_mat_without_cp_1C, N_SC, 1);
syms_f_mat_1D = fft(rx_payload_mat_without_cp_1D, N_SC, 1);
    
                                                                %% MRC combining of 1x4 RX symbols
%syms_f_eq_1A = syms_f_mat_1A ./ transpose(H_cap_A1A_f);
% syms_f_eq_1A = syms_f_mat_1A ./ transpose(H_cap_A1A_f);
% syms_f_eq_1A = syms_f_mat_1A ./ transpose(H_cap_A1A_f);
% syms_f_eq_1A = syms_f_mat_1A ./ transpose(H_cap_A1A_f);

rx_1x4_mrc = syms_f_mat_1A .* conj(transpose(H_cap_A1A_f)) + syms_f_mat_1B .* conj(transpose(H_cap_A1B_f)) ...
                            + syms_f_mat_1C .* conj(transpose(H_cap_A1C_f)) + syms_f_mat_1D .* conj(transpose(H_cap_A1D_f)) ;
rx_1x4_mrc = rx_1x4_mrc ./(transpose(abs(H_cap_A1A_f).^2) + ... 
                                                          transpose(abs(H_cap_A1B_f).^2) + ...
                                                          transpose(abs(H_cap_A1C_f).^2) + ...
                                                          transpose(abs(H_cap_A1D_f).^2));


                                                                      %% 1x4 MIMO SFO correction 

if DO_APPLY_SFO_CORRECTION
    sfo_coorected_rx_Sig =sfo(rx_1x4_mrc,N_OFDM_SYMS,N_SC);
else 
    sfo_coorected_rx_Sig = rx_1x4_mrc;
end
    
                                                                    %% 1x4 MIMO Phase correction 

if DO_APPLY_PHASE_ERR_CORRECTION
    phase_corrected_1x4_signal = phase_correction(sfo_coorected_rx_Sig,pilots_A,N_OFDM_SYMS,N_SC);
else 
    phase_corrected_1x4_signal = sfo_coorected_rx_Sig;
end
                                                                 %% 1x4 MIMO Demodulate and Demap

post_processed_1x4_mimo_rx_sig_data_only = phase_corrected_1x4_signal(SC_IND_DATA,:);
reshaped_rx_signal_1x4 = reshape(post_processed_1x4_mimo_rx_sig_data_only,1,(N_SC-16)*N_OFDM_SYMS);

% figure;
% scatter(real(reshaped_rx_signal_1x4), imag(reshaped_rx_signal_1x4),'filled');
% title(' 1x4 MIMO Combining Signal Space of received bits');
% xlabel('I'); ylabel('Q');

demap_1x4_mimo = demapper(reshaped_rx_signal_1x4, MOD_ORDER, 1);
decoded_rx_signal_1x4= vitdec(demap_1x4_mimo,trel,7,'trunc','hard');
[number_1x4,ber_1x4] = biterr(tx_data_a,decoded_rx_signal_1x4);
BER1x4 = ber_1x4;

end


if mimo_2x2
            
                                                                                              % 2x2 MIMO Payload extraction %
        
        rx_payload_mat_with_cp_2A = reshape(rx_dec_cfo_corr_2A(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);
        rx_payload_mat_with_cp_2B = reshape(rx_dec_cfo_corr_2B(1,payload_ind+1:payload_ind+(N_SC+CP_LEN)*N_OFDM_SYMS), N_SC+CP_LEN, N_OFDM_SYMS);

                                                                                                           % CP removal %
        rx_payload_mat_without_cp_2A = rx_payload_mat_with_cp_2A(CP_LEN+[1:N_SC]:end,:);
        rx_payload_mat_without_cp_2B = rx_payload_mat_with_cp_2B(CP_LEN+[1:N_SC]:end,:);

                                                                                                         % FFT 1x4 MIMO %
        syms_f_mat_2A = fft(rx_payload_mat_without_cp_2A);
        syms_f_mat_2B = fft(rx_payload_mat_without_cp_2B);
    
                                                                                                %% 2x2 MIMO Combining %%

         rx_2x2_mimo_combining = (conj(transpose(H_cap_11_A_f+H_cap_12_A_f )) .* syms_f_mat_2A) + ...
                                                            (conj(transpose(H_cap_22_B_f + H_cap_21_B_f)) .* syms_f_mat_2B);

        rx_2x2_mimo_combining = rx_2x2_mimo_combining ./ (transpose(abs(H_cap_11_A_f+H_cap_12_A_f).^2)...
                                                                                                                + transpose(abs(H_cap_22_B_f + H_cap_21_B_f).^2)) ;

                                                                                                 %% 2x2 MIMO SFO correction 

    if DO_APPLY_SFO_CORRECTION
        sfo_coorected_rx_Sig_2x2 =sfo(rx_2x2_mimo_combining,N_OFDM_SYMS,N_SC);
    else 
        sfo_coorected_rx_Sig_2x2 = rx_2x2_mimo_combining;
    end

                                                                                                 %% 2x2  MIMO Phase correction 
    if DO_APPLY_PHASE_ERR_CORRECTION
        phase_corrected_2x2_signal = phase_correction(sfo_coorected_rx_Sig_2x2,pilots_A,N_OFDM_SYMS,N_SC);
    else 
        phase_corrected_2x2_signal = sfo_coorected_rx_Sig_2x2;
    end

                                                                                            %% 2x2 MIMO Demodulate and Demap

    post_processed_2x2_mimo_rx_sig_data_only = phase_corrected_2x2_signal(SC_IND_DATA,:);
    reshaped_rx_signal_2x2 = reshape(post_processed_2x2_mimo_rx_sig_data_only,1,(N_SC-16)*N_OFDM_SYMS);
    
%     figure;
%     scatter(real(reshaped_rx_signal_2x2), imag(reshaped_rx_signal_2x2),'filled');
%     title(' 2x2 MIMO Combining Signal Space of received bits');
%     xlabel('I'); ylabel('Q');
    
    demap_2x2_mimo = demapper(reshaped_rx_signal_2x2, MOD_ORDER, 1);
    decoded_rx_signal_2x2= vitdec(demap_2x2_mimo,trel,7,'trunc','hard');
    [number_2x2,ber_2x2] = biterr(tx_data_a, decoded_rx_signal_2x2);
    BER2x2 = ber_2x2;
end

end
