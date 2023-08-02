function errors=qam_errors(k,Nsymb,nsamp,EbNo)
% Εξομοιώνεται το πλήρες σύστημα πομπού-δέκτη QAM, με αφετηρία
% δυαδική ακολουθία προς εκπομπή. Γίνεται κωδικοποίηση
% Gray, ενώ Γκαουσιανός θόρυβος προστίθεται στο ζωνοπερατό
% QAM σήμα για δοσμένο σηματοθορυβικό λόγο Eb/No.
% Δίνεται δυνατότητα επιλογής φίλτρου μορφοποίησης βασικής ζώνης
% (ορθογωνικό ή Nyquist), ενώ μετράται το BER μετά τη βέλτιστη
% αποδιαμόρφωση με το αντίστοιχο προσαρμοσμένο φίλτρο.
% M (=2^k, k άρτιο)είναι το μέγεθος του σηματικού αστερισμού,
% mapping είναι το μιγαδικό διάνυσμα των Μ QAM συμβόλων,
% σε διάταξη κωδικοποίησης Gray
% Nsymb είναι το μήκος της ακολουθίας QAM,
% nsamp είναι ο συντελεστής υπερδειγμάτισης, (# δειγμάτων/Τ)
% EbNo είναι ο λόγος Eb/No, in db
%
% LxL σηματικός αστερισμός, L=2^l, l>0
M=2^k;
L=sqrt(M);
l=log2(L);
pulse_type=1; % επιλογή rtrapezium ως φίλτρου μορφοποίησης
% pulse_type=0; % επιλογή μορφοποίησης ορθογωνικού παλμού
fc=10; % συχνότητα φέροντος, πολλαπλάσιο του Baud Rate (1/T)
SNR=EbNo-10*log10(nsamp/k/2); % SNR ανά δείγμα σήματος
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
core=[1+i;1-i;-1+i;-1-i];
mapping=core;
if(l>1)
for j=1:l-1
mapping=mapping+j*2*core(1);
mapping=[mapping;conj(mapping)];
mapping=[mapping;-conj(mapping)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ΠΟΜΠΟΣ %%%%%%%%%%%%
x=floor(2*rand(k*Nsymb,1)); % τυχαία δυαδική ακολουθία
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb')';
y=[];
for n=1:length(xsym)
y=[y mapping(xsym(n)+1)];
end
%% Ορισμός φίλτρου μορφοποίησης
if (pulse_type==1) % παλμός Nyquist -- rtrapezium
delay = 8; % Group delay (# περιόδων Τ)
filtorder = delay*nsamp*2;
rolloff = 0.25; % συντελεστής εξάπλωσης φίλτρου
% Η rtrapezium.m πρέπει να βρίσκεται στο current directory
shaping_filter = rtrapezium(nsamp,rolloff,delay);
else % ορθογωνικός παλμός
delay=0.5;
shaping_filter=ones(1,nsamp)/sqrt(nsamp);% με κανονικοποίηση!
end
%% Εκπεμπόμενο σήμα
ytx=upsample(y,nsamp);
ytx = conv(ytx,shaping_filter);
% Υπολογισμός & σχεδιασμός φάσματος
%figure(1); pwelch(real(ytx),[],[],[],nsamp); % σε κλίμακα 1/T
% quadrature modulation
m=(1:length(ytx));
s=real(ytx.*exp(1j*2*pi*fc*m/nsamp));
%figure(2); pwelch(s,[],[],[],nsamp); % Uncomment for lab5.(2) σε κλίμακα 1/T
% ---------------------
% Διάγραμμα οφθαλμού, όταν ορθ. παλμός.
% if (pulse_type==1) eyediagram(ytx(1:2000),nsamp*2);
% Προσθήκη λευκού γκαουσιανού θορύβου
Ps=10*log10(s*s'/length(s)); % ισχύς σήματος, σε db
Pn=Ps-SNR; % αντίστοιχη ισχύς θορύβου, σε db
n=sqrt(10^(Pn/10))*randn(1,length(ytx));
snoisy=s+n; % θορυβώδες ζωνοπερατό σήμα
clear ytx xsym s n; % για εξοικονόμηση μνήμης
%% ΔΕΚΤΗΣ %%
% Αποδιαμόρφωση
yrx=2*snoisy.*exp(-1j*2*pi*fc*m/nsamp); clear s;
% Κανονικά ακολουθεί βαθυπερατό φίλτρο.
% Όμως στην περίπτωση μορφοποίησης Nyquist δεν χρειάζεται,
% αφού το προσαρμοσμένο φίλτρο (Nyquist) είναι ταυτόχρονα
% και ένα καλό βαθυπερατό φίλτρο.
yrx = conv(yrx,shaping_filter);
%figure(3); pwelch(real(yrx),[],[],[],nsamp); % κλίμακα 1/T
yrx = downsample(yrx,nsamp); % Υποδειγμάτιση στο πλέγμα nT.
yrx = yrx(2*delay+(1:length(y))); % περικοπή άκρων συνέλιξης.
% ----------------------
yi=real(yrx); yq=imag(yrx); % συμφασική και εγκάρσια συνιστώσα
xrx=[]; % διάνυσμα δυαδικής ακολουθίας εξόδου -- αρχικά κενό
q=[-L+1:2:L-1];
for n=1:length(yrx) % επιλογή πλησιέστερου σημείου
[m,j]=min(abs(q-yi(n)));
yi(n)=q(j);
[m,j]=min(abs(q-yq(n)));
yq(n)=q(j);
end
err=not(y==(yi+i*yq));
errors=sum(err);
%scatterplot(yrx); % διάγραμμα διασκορπισμού
% Τι παρατηρείτε στο διάγραμμα διασκορπισμού, στην περίπτωση
% μορφοποίησης με ορθ. παλμό; Ποια η επίπτωση στο ber
% και πώς διορθώνεται;
%% ΕΚΤΙΜΗΣΗ ΡΥΘΜΟΥ ΛΑΘΩΝ
% Με δύο τρόπους: (α) με σύγκριση των σημείων QAM (y-yrx)
% (β) με σύγκριση των δυαδικών ακολουθιών (x-xrx)
%gber1=sum(not(y==(yi+i*yq)));
%ber2=sum(not(xrx==x));
%errors=ber2;
end

