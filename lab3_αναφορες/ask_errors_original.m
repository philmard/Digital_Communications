function errors=ask_errors_original(k,M,nsamp,EbNo)
2. % Η συνάρτηση αυτή εξομοιώνει την παραγωγή και αποκωδικοποίηση
3. % θορυβώδους σήματος L-ASK και μετρά τον αριθμό των εσφαλμένων συμβόλων.
% Υπολογίζει επίσης τη θεωρητική πιθανότητα εσφαλμένου συμβόλου, Pe.
% Επιστρέφει τον αριθμό των εσφαλμένων συμβόλων, καθώς και τον συνολικό
% αριθμό των συμβόλων που παρήχθησαν.
% k είναι ο αριθμός των bits/σύμβολο, ώστε L=2^k,
% M ο αριθμός των παραγόμενων συμβόλων (μήκος ακολουθίας L-ASK)
% nsamp ο αριθμός των δειγμάτων ανά σύμβολο (oversampling ratio)
% EbNo είναι ο λόγος Eb/No, σε db
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
% Διάνυσμα τυχαίων ακεραίων {±1, ±3, ... ±(L-1)}. Να επαληθευτεί
x=2*floor(L*rand(1,M))-L+1;
Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
Pmetr=sum(x.^2)/length(x); % μετρούμενη ισχύς σήματος (για επαλήθευση)
y=rectpulse(x,nsamp);
n=wgn(1,length(y),10*log10(Px)-SNR);
ynoisy=y+n; % θορυβώδες σήμα
y=reshape(ynoisy,nsamp,length(ynoisy)/nsamp);
matched=ones(1,nsamp);
z=matched*y/nsamp;
l=[-L+1:2:L-1];
for i=1:length(z)
[m,j]=min(abs(l-z(i)));
z(i)=l(j);
end
err=not(x==z);
errors=sum(err);
end