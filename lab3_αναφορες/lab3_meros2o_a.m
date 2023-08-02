L = 16;
EbNodB = 0:0.1:20;
EbNo = 10.^(EbNodB/10);
Pe = ((L-1)/L)*erfc(sqrt(3*log2(L)*EbNo/((L^2)-1)));
BER = Pe/log2(L);
semilogy(EbNodB, BER, '--r', 'LineWidth', 2)
xlabel('Eb/No (dB)')
ylabel('BER')

Nerrors10 = ask_errors_original(4,20000,20,10); % L=16
Nerrors15 = ask_errors_original(4,20000,20,15); % L=16
Nerrors20 = ask_errors_original(4,20000,20,20); % L=16
M=20000;
Pe10 = Nerrors10/M;
Pe15 = Nerrors15/M;
Pe20 = Nerrors20/M;
BER10 = Pe10/log2(L)
BER15 = Pe15/log2(L)
BER20 = Pe20/log2(L)

