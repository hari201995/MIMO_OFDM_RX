function sfo_corrected_signal = sfo(rx_sig, N_OFDM_SYMS, N_SC)
    sfo_corrected_signal = zeros(N_SC,N_OFDM_SYMS);
    % extract pilots 
    SC_IND_PILOTS           = [8 22 44 58];  
    rx_pilot = rx_sig(SC_IND_PILOTS,:);
    pilot_phase_change_across_sym = zeros(4,N_OFDM_SYMS);
    % pilot phase difference across symbols 
    for i=1:N_OFDM_SYMS-1
    
        delta_phi_1 = angle(rx_pilot(1,i+1)) - angle(rx_pilot(1,i));
        delta_phi_2 = angle(rx_pilot(2,i+1)) - angle(rx_pilot(2,i));
        delta_phi_3 = angle(rx_pilot(3,i+1)) - angle(rx_pilot(3,i));
        delta_phi_4 = angle(rx_pilot(4,i+1)) - angle(rx_pilot(4,i));
       
        del_phi = [delta_phi_1; delta_phi_2; delta_phi_3; delta_phi_4];
        pilot_phase_change_across_sym(:,i+1) = del_phi;
    end
    
    sfo_corrected_signal(:,1) = rx_sig(:,1);
    for i=2:N_OFDM_SYMS
        sfo_corrected_signal(:,i) =  rx_sig(:,i) .* exp(-1i*mean(pilot_phase_change_across_sym(:,i)));
    end
end