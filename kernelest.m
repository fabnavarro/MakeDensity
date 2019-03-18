function [y,rot] = kernelest(x,t,kname,h)
%  Standard kernel density estimation
%  (the rule-of-thumb is used by default if a bandwidth is not specified)
% 	Input
%		x     Input realizations.
%       t     Vector ti of values where the density estimate is evaluated.
%       h     The bandwidth of the kernel smoothing window.
%           
%	Output
%		y     Estimated Density.
%       rot   The ROT bandwidth of the kernel smoothing window used.
%             (The default is optimal for estimating normal densities) 
if (exist('x')~=1), error('Provide input realizations'); end;
if (exist('t')~=1), error('Provide t'); end;
if (exist('kname')~=1 || isempty(kname)), kname = 'Normal'; end;
m = length(x);
if (exist('h')~=1 || isempty(h)), h = rotbandw(x,m,kname); end;

n = length(t);
K = @(x)kerneltype(x,kname);
y = zeros(1,n);
for j = 1:n
  y(j) = sum(K((t(j)-x)/h));
end
y = y/(m*h);
y = y(:);
rot = rotbandw(x,m,kname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = kerneltype(x,Name)
%  Some kernels
% 	Input
%		x Input vector.
%       Names:  normal,
%       p = 0: Uniform kernel.
%       p = 1: Epanechnikov kernel.
%       p = 2: Biweight kernel.
B = @(a,b)(gamma(a)*gamma(b)/gamma(a+b));
if strcmp(Name,'Normal'),
  s = 1; 
  m = 0;
  k = exp(-0.5*((x-m)./s).^2)./(sqrt(2*pi)*s); 
elseif strcmp(Name,'Uniform'),
  p = 0;  
  k = (1-x.^2).^p/(2^(2*p+1)*B(p+1,p+1)).*(abs(x)<1);
elseif strcmp(Name,'Epanechnikov'),
  p = 1;
  k = (1-x.^2).^p/(2^(2*p+1)*B(p+1,p+1)).*(abs(x)<1);
elseif strcmp(Name,'Biweight'),
  p = 2;
  k = (1-x.^2).^p/(2^(2*p+1)*B(p+1,p+1)).*(abs(x)<1);  
else
  disp('Unkwnown kernel type.');
  disp('Allowable Names are:')
  disp('Normal'),
  disp('Uniform'),
  disp('Epanechnikov'),
  disp('Biweight'),
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rot = rotbandw(x,n,Name)
med = median(x);
sig = median(abs(x-med))/0.6745;
if strcmp(Name,'Normal'),
  rot = sig*(4/(3*n))^(1/5);
elseif strcmp(Name,'Uniform'),
  rot = sig*(4/(3*n))^(1/5);
elseif strcmp(Name,'Epanechnikov'),
  rot = 2.34*sig*n^(-1/5);
elseif strcmp(Name,'Biweight'),        
  rot = 2.78*sig*n^(-1/5);
else
  disp('Unkwnown kernel type.');
  disp('Allowable Names are:')
  disp('Normal'),
  disp('Uniform'),
  disp('Epanechnikov'),
  disp('Biweight'),
end