T = 1000;
n = 5000;
name = 'Claw';
[pdfx,x,t] = MakeDensity(name,T,n);
[y,~] = kernelest(x,t);
figure
plot(t,pdfx,t,y,'linewidth',2);axis tight
title(name)
