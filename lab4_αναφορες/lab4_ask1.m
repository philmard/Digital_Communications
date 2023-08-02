clear all; close all; clc;
k=3; Nsymb=3333; nsamp=32;
L=2^k;
x=round(rand(1,k*Nsymb)); % παραγωγη τυχαιας δυαδικης ακολουθιας
step=2; % ο αριθμός των πλατών και το βήμα μεταξύ τους k=log(L); 
mapping=[step/2; -step/2]; 
if(k>1) 
    for j=2:k 
        mapping=[mapping+2^(j-1)*step/2; ... 
            -mapping-2^(j-1)*step/2]; 
    end 
end 
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb'); 
x=[]; 
for i=1:length(xsym) 
    x=[x mapping(xsym(i)+1)]; 
end
% Ορισμός παραμέτρων φίλτρου
nsamp=32;
delay = 8; % Group delay (# of input symbols)>5 για καλύτερα αποτελέσματα
filtorder = delay*nsamp*2; % τάξη φίλτρου
rolloff = 0.4; % Συντελεστής πτώσης -- rolloff factor 0.1 ή 0.9
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
y=upsample(x,nsamp);
ytx = conv(y,rNyquist); 
%% Λαμβανόμενο σήμα: ytx (χωρίς παραμόρφωση)
% Φιλτράρισμα σήματος με φίλτρο τετρ. ρίζας ανυψ. συνημ.
yrxf=conv(ytx,rNyquist);
% ----------------------
yrx = downsample(yrxf,nsamp);
yrx=yrxf(2*delay*nsamp+1:end-2*delay*nsamp); %Περικοπή-λόγω καθυστέρησης
% Σχεδίαση yrx και υπέρθεση x?
figure(1)
plot(yrx(1:10*nsamp)); %Για να σχεδιάσουμε το σήμα σε 10Τ
hold;
stem(y(1:10*nsamp),'filled');%Για να υπερθέσουμε σε Τ=nsamp
grid;
figure(2)
pwelch(yrx);