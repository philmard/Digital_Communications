clear all; close all;
% Το αρχείο "sima.mat" περιέχει το σήμα s και τη συχνότητα
% δειγματοληψίας Fs. Το φάσμα του σήματος εκτείνεται σχεδόν σε όλη την
% περιοχή συχνοτήτων μέχρι 4 KHz. Πάνω από 1 KHz, όμως, είναι θόρυβος
% και πρέπει να φιλτραριστεί.
load sima;
figure; pwelch(s,[],[],[],Fs);
f1=700; f2=1500;
Ts=1/Fs;
f2m1=(f2-f1); f2p1=(f2+f1)/2; N=256; %
t=[-(N-1):2:N-1]*Ts/2; %#ok<*NBRAK1>
hbp=2/Fs*cos(2*pi*f2p1*t).*sin(pi*f2m1*t)/pi./t;
hbpw=hbp.*kaiser(length(hbp),5)';
sima_bp_k=conv(s,hbpw);
figure; pwelch(sima_bp_k,[],[],[],Fs);
%hbppm = firpm(64, [0 0.10 0.15 0.5]*2, [1 1 0 0]);
f=2*[0 f1*0.95 f1*1.05 f2*0.95 f2*1.05 Fs/2]/Fs;
hbp_pm=firpm(256, f, [0 0 1 1 0 0]);
% figure; freqz(hpm,1);
wvtool(hbp_pm);
sima_bp_pm=conv(s,hbp_pm);
figure; pwelch(sima_bp_pm,[],[],[],Fs);
sound(20*s); % ακούμε το αρχικό σήμα, s
pause;
sound(20*sima_bp_pm); % ακούμε το φιλτραρισμένο σήμα, sima_bp me parks
pause;
sound(20*sima_bp_k); % ακούμε το φιλτραρισμένο σήμα, sima_bp me kaiser