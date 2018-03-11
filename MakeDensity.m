function [pxtheo,x,t] = MakeDensity(Name,T,n)
%
%  Make artificial normal mixture density and
%  generate associated normal mixture distributed random numbers for
%  the fifteen density examples used Marron and Wand, AOS, 1992.
%
% 	Input
%		Name   String: 'Gaussian','SkewedUnimodal','StronglySkewed',
%                      'Kurtotic','Outlier','Bimodal','SeparatedBimodal',
%                      'AsymmetricBimodal','Trimodal','Claw','DoubleClaw',
%                      'AsymmetricClaw','AsymmetricDoubleClaw'
%                      'SmoothComb','DiscreteComb'
%                       (Wand & Marron densities)
%		n	   Desired sampling length.
%		T	   Output density length.
%
%	Output
%		pxtheo Theoretical Density.
%       x      Random samples.
%  
%  Copyright (c) 2018 Fabien Navarro

if nargin > 1,
  t = linspace(-3,3,T)';
end
if strcmp(Name,'Gaussian'),
  q = 1;
  m = 0;
  s = 1;
elseif  strcmp(Name,'SkewedUnimodal'),
  q = [1/5 1/5   3/5];
  m = [0   1/2 13/12];
  s = [1   2/3   5/9];
elseif  strcmp(Name,'StronglySkewed'),
  q = [1   1    1    1    1      1      1       1       ]/8;
  m = [0  -1   -5/3 -19/9 -65/27 -633/243 -1995/729 -6177/2187 ];
  s = [1   2/3  4/9  8/27 16/81  32/243 64/729  128/2187];
elseif  strcmp(Name,'Kurtotic'),
  q = [ 2 1]/3;
  m = [ 0 0];
  s = [10 1]/10;
elseif  strcmp(Name,'Outlier'),
  q = [ 1 9]/10;
  m = [ 0 0];
  s = [10 1]/10;
elseif  strcmp(Name,'Bimodal'),
  q = [ 1 1]/2;
  m = [-1 1];
  s = [ 2 2]/3;
elseif  strcmp(Name,'SeparatedBimodal'),
  q = [ 1 1]/2;
  m = [-3 3]/2;
  s = [ 1 1]/2;
elseif  strcmp(Name,'AsymmetricBimodal'),
  q = [3 1]/4;
  m = [0 3/2];
  s = [1 1/3];
elseif  strcmp(Name,'Trimodal'),
  q = [ 9/20 9/20 1/10];
  m = [-6/5  6/5  0   ];
  s = [ 3/5  3/5  1/4 ];
elseif  strcmp(Name,'Claw'),
  q = [1/2 1/10 1/10 1/10 1/10 1/10];
  m = [0 -1 1/2-1 0 3/2-1 1 5/2-1];
  s = [1 1/10 1/10 1/10 1/10 1/10];
elseif  strcmp(Name,'DoubleClaw'),
  q = [ 49/100 49/100 1/350 1/350 1/350 1/350 1/350 1/350 1/350];
  m = [-1      1     -3/2  -1    -1/2   0     1/2   1     3/2  ];
  s = [ 2/3    2/3    1/100 1/100 1/100 1/100 1/100 1/100 1/100];
elseif  strcmp(Name,'AsymmetricClaw'),
  q = [ 1/2  8/31  4/31 2/31 1/31 1/62];
  m = [ 0   -3/2  -1/2  1/2  3/2  5/2];
  s = [ 1    4/10  2/10 1/10 1/20 1/40];
elseif  strcmp(Name,'AsymmetricDoubleClaw'),
  q = [  46/100 46/100  1/300 1/300 1/300 7/300 7/300 7/300];
  m = [ -1      1      -1/2  -1    -3/2   1/2   1     3/2  ];
  s = [  2/3    2/3     1/100 1/100 1/100 7/100 7/100 7/100];
elseif  strcmp(Name,'SmoothComb'),
  q = [ 32 16  8  4  2  1]/63;
  m = [-31 17 41 53 59 62]/21;
  s = [ 32 16  8  4  2  1]/63;
elseif  strcmp(Name,'DiscreteComb'),
  q = [ 2/7   2/7 2/7 1/21  1/21  1/21];
  m = [-15/7 -3/7 9/7 16/7 18/7 20/7];
  s = [ 2/7   2/7 2/7 1/21  1/21  1/21];
else
  disp('Unkwnown Density type.');
  disp('Allowable Names are:')
  disp('Gaussian'),
  disp('SkewedUnimodal'),
  disp('StronglySkewed'),
  disp('Kurtotic'),
  disp('Outlier'),
  disp('AsymmetricBimodal'),
  disp('Bimodal'),
  disp('SeparatedBimodal'),
  disp('Trimodal'),
  disp('Claw'),
  disp('DoubleClaw'),
  disp('AsymmetricClaw'),
  disp('AsymmetricDoubleClaw'),
  disp('SmoothComb'),
  disp('DiscreteComb')
end
pdfx   = @(x,q,m,s)dgauss(x, q, m, s);
pxtheo = pdfx(t',q,m,s)';
x      = MixGauss(n,q,m,s)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  f = dnorm(x,m,s)
%DNORM 	  The normal density function
%
%         f = dnorm(x,Mean,StandardDeviation)

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg

if nargin<3, s=1; end
if nargin<2, m=0; end
f = exp(-0.5*((x-m)./s).^2)./(sqrt(2*pi)*s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  X = rnorm(n,m,s)
%RNORM 	  Normal random numbers
%
%         p = rnorm(Number,Mean,StandardDeviation)

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg

if nargin<3, s=1; end
if nargin<2, m=0; end
if nargin<1, n=1; end
if size(n)==1
  n = [n 1];
  if size(m,2)>1, m = m'; end
  if size(s,2)>1, s = s'; end
end
X = randn(n).*s + m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = dgauss(D,p,m,s)
% gaussian mixture density
d = zeros(1,length(D));
for j = 1:length(p)
  d = d + p(j)*dnorm(D,m(j),s(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = MixGauss(n,p,m,s)
% iid sample of random variables from gaussian mixture density
y = rand(1,n);
X = zeros(1,n);
q = [0 cumsum(p)];
for j = 1:length(p)
  ind = (y<=q(j+1)) & (y>q(j));
  X(ind) = rnorm(sum(ind+0), m(j), s(j));
end