function [b,Cb,vf]=lreg(x,y,Cy)
% Linear Regression
%   Computes slope and y-intercept for a linear
%   regression of y on x.  Can optionally weight
%   observations y.
% Version 17 Jun 97
% Usage:  [b,Cb,vf]=lreg(x,y)    -> simple regression
%         [b,Cb,vf]=lreg(x,y,Cy) -> weighted regression
% Input:  x  - vector of independent variable
%         y  - vector of dependent variable
%         Cy - a priori covariance matrix for y
%              (default = identity matrix)
% Output: b  - estimated parameters
%              b(1) = intercept
%              b(2) = slope
%         Cb - covariance matrix for b (scaled by vf)
%         vf - estimated variance factor
if nargin==1
  error('Incorrect number of arguments');
end
n=max(size(y));
if nargin==2
  A=[ones(n,1) x];
  Cb=inv(A'*A);
  b=Cb*A'*y;
  r=y-A*b;
  vf=(r'*r)/(n-2);
elseif nargin==3
  P=inv(Cy);
  A=[ones(n,1) x];
  AtP=A'*P;
  Cb=inv(AtP*A);
  b=Cb*AtP*y;
  r=y-A*b;
  vf=(r'*P*r)/(n-2);
end
Cb=vf*Cb;
