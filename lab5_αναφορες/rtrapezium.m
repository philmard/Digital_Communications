% function rtrapezium -- square root trapezium (Nyquist filter)
%
function [rtrfilter Ho]=rtrapezium(nsamp,rolloff,delay)
F0 = 1/2/nsamp; %
% Υπερδειγμάτιση, για λόγους ακρίβειας
dense=32;
nsamp1=nsamp*dense;
a = rolloff;
F1=F0*(1-a); F2=F0*(1+a); % ακραίες συχνότητες
N = 2*delay*nsamp1; % τάξη φίλτρου, εδώ μήκος κρουστικής απόκρισης
for k=1:N/2
f=(k-1)/N;
if (f<F1) Ho(k)=1;
elseif (f>F2) Ho(k)=0;
else Ho(k)=max(0,1-(f-F1)/(F2-F1));
end
Ho(N/2+1)=0;
Ho(N+2-k)=Ho(k);
end
H=sqrt(Ho).*exp(j*pi*(N+2)/(N+1)*[0:N]);
h=ifft(H,'symmetric');
% περικοπή ουράς (ορθ. παράθυρο)
M=N/dense;
for k=-M/2:M/2
rtrfilter(k+M/2+1)=h(N/2+1+k);
end
% κανονικοποίηση
rtrfilter=rtrfilter/sqrt(sum(rtrfilter.^2));
end