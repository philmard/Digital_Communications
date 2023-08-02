function errors=ask_Nyq_filter(k,Nsymb,nsamp,EbNo)
% Η συνάρτηση αυτή εξομοιώνει την παραγωγή και αποκωδικοποίηση
% θορυβώδους ακολουθίας L-ASK και μετρά τα λαθεμένα σύμβολα,
% με μορφοποίηση παλμών μέσω φίλτρου τετρ. ρίζας Nyquist.
%%%%%% ΠΑΡΑΜΕΤΡΟΙ %%%%%%%
% k είναι ο αριθμός των bits ανά σύμβολο, έτσι L=2^k
% Nsymb είναι το μήκος της εξομοιούμενης ακολουθίας συμβόλων L-ASK
% nsamp είναι ο συντελεστής υπερδειγμάτισης, δηλ. #samples/Td
% EbNo είναι ο ανηγμένος σηματοθορυβικός λόγος, Εb/No, σε db
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
% Κωδικοποίηση Gray σύμφωνα με την εκφώνηση
step=2; %το βήμα μεταξύ των πλατών
x=randi(2,1,k*Nsymb)-1; 
mapping=[step/2; -step/2]; 
if(k>1) 
   for j=2:k 
         mapping=[mapping+2^(j-1)*step/2; ... % gray
                 -mapping-2^(j-1)*step/2];    % gray     
         %mapping=-(L-1):step:(L-1); % NOT gray
   end 
end 
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb'); 
y=[];
for i=1:length(xsym) 
    y=[y mapping(xsym(i)+1)]; 
end
%% Ορισμός παραμέτρων φίλτρου 
%delay = 1;% Group delay 
%delay=2;
delay=5;
filtorder = delay*nsamp*2; % τάξη φίλτρου
rolloff = 0.1;% rolloff factor
%rolloff=0.2;
%rolloff=0.4;
%rolloff=0.5; %Για το 4ο ερώτημα

% κρουστική απόκριση φίλτρου τετρ. ρίζας ανυψ. συνημιτόνου
rNyquist= rcosine(1,nsamp,'fir/sqrt',rolloff,delay);
% ----------------------
% Για φίλτρο γραμμικής πτώσης να χρησιμοποιηθεί η επόμενη εντολή
% (με την rtrapezium του Κώδικα 4.1 στο current directory)
% Πρέπει delay>5 για καλά αποτελέσματα στην αναγνώριση.
% rNyquist=rtrapezium(nsamp,rolloff,delay);
% ----------------------
%% ΕΚΠΕΜΠΟΜΕΝΟ ΣΗΜΑ
% Υπερδειγμάτιση και εφαρμογή φίλτρου rNyquist
y=upsample(y,nsamp);
ytx = conv(y,rNyquist);
ynoisy=awgn(ytx,SNR,'measured'); % θορυβώδες σήμα
% ----------------------
%% ΛΑΜΒΑΝΟΜΕΝΟ ΣΗΜΑ
% Φιλτράρισμα σήματος με φίλτρο τετρ. ρίζας ανυψ. συνημ.
yrx=conv(ynoisy,rNyquist);
yrx = downsample(yrx,nsamp); % Υποδειγμάτιση
yrx = yrx(2*delay+1:end-2*delay); % περικοπή, λόγω καθυστέρησης
% Ανιχνευτής ελάχιστης απόστασης L πλατών
%l=[-L+1:step:L-1];
Xr=[];
for i=1:length(yrx)
[m,j]=min(abs(mapping-yrx(i)));
xr=de2bi(j-1,k,'left-msb');
Xr=[Xr xr];
end
err=not(x==Xr);
errors=sum(err);


