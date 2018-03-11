T = 1000;
n = 5000;
name = 'Gaussian';
[pdfx,~,t] = MakeDensity(name,T,n);       
figure
plot(t,pdfx,'linewidth',2);axis tight
title(name)