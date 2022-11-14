function phase_corrected_signal = phase_correction(rx_sig, pilots_A,N_OFDM_SYMS, ~)
    SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
    %extract pilots 
    SC_IND_PILOTS           = [8 22 44 58];  
    rx_pilot = rx_sig(SC_IND_PILOTS,:);
    phase_corrected_signal = rx_sig;
    for i=1:N_OFDM_SYMS
        phase_diff_pilot_1 = angle (rx_pilot(1,i) / pilots_A(1));
        phase_diff_pilot_2 = angle (rx_pilot(2,i) / pilots_A(2));
        phase_diff_pilot_3 = angle (rx_pilot(3,i) / pilots_A(3));
        phase_diff_pilot_4 = angle (rx_pilot(4,i) / pilots_A(4));
        
        phase_correction = mean([phase_diff_pilot_1 phase_diff_pilot_2 phase_diff_pilot_3 phase_diff_pilot_4 ]);
        phase_corrected_signal(SC_IND_DATA,i) = rx_sig(SC_IND_DATA,i) *exp(-1i*phase_correction);
    
    end
end