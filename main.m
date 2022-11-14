close all
clear all

ITR=50; % number of iteration
snr = 1:1:20; % snr values.
error_mat_1x4 = zeros(ITR,length(snr));
error_mat_2x2= zeros(ITR,length(snr));

for idx=1:length(snr)
    for itr = 1:50
        [ber1x4, ber2x2]= MIMO_OFDM_RX_test(snr(idx));
        error_mat_1x4(itr,idx) = ber1x4;
        error_mat_2x2(itr,idx) = ber2x2;
        X = sprintf('%d ITR - %d SNR \n',itr, snr(idx));
        disp(X);
    end
end


for i=1:length(snr)
variance_ber_1x4(i) = var(error_mat_1x4(:,i));
mean_ber_1x4(i) = mean(error_mat_1x4(:,i));
variance_ber_2x2(i) = var(error_mat_2x2(:,i));
mean_ber_2x2(i) = mean(error_mat_2x2(:,i));
end

plot(snr,variance_ber_1x4);
hold on;
plot(snr,variance_ber_2x2);
plot(mean_ber_1x4);
plot(mean_ber_2x2);

 legend('ber variance 1x4 mimo', 'ber variance 2x2 mimo','ber mean 1x4 mimo','ber mean 2x2 mimo');
title('BER mean and variance across snr');
xlabel('snr');
ylabel('BER');

