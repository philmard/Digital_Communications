close all;
clear all;

M=64;
L=8;
l=3;

core=[1+i;1-i;-1+i;-1-i]; % τετριμμένη κωδικοποίηση, M=4
mapping=core;
if(l>1)
 for j=1:l-1
 mapping=mapping+j*2*core(1);
 mapping=[mapping;conj(mapping)];
 mapping=[mapping;-conj(mapping)];
 end
end

% Plot the constellation points
scatterplot(mapping);
hold on;

% Add binary words near each point
for i=1:length(mapping)
text(real(mapping(i)),imag(mapping(i)),num2str(de2bi(i-1,2*l,'left-msb')),'FontSize', 6);
end

hold off;

