L = 4;
EbNo = 0:0.1:20;
Pe = ((L-1)/L)*erfc(sqrt(3*log2(L)*EbNo/((L^2)-1)));
k = log2(L);
BER = Pe/(k * length(EbNo));
semilogy(EbNo, BER, '--r', 'LineWidth', 2)
xlabel('Eb/No (dB)')
ylabel('BER')