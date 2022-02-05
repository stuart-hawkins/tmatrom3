% Gauss-Legendre quadrature rule
%
%   [p,w]=gauss_legendre(n) returns points p and weights w for the
%   n+1 point Gauss-Legendre rule.
%
% See also: legendre, gauss_hermite.
%
% Stuart C. Hawkins - 19 February 2015

% Copyright 2019-2022 Stuart C. Hawkins
% 	
% This file is part of TMATROM3
% 
% TMATROM3 is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM3 is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM3.  If not, see <http://www.gnu.org/licenses/>.


function [pp,pw]=gauss_legendre(n)

% upper limit for the seach interval for the roots
my_inf=1;

% set root of first polynomial
% Note: roots{k} gives roots of the degree k-1 polynomial.
roots{1}=[ 0 ];

% loop through polynomial degree
for k=1:n
    
    % we use the interlace property to find the roots of the degree k
    % polynomial from the roots of the degree k-1 polynomial    
    ints=[-my_inf;roots{k};my_inf];
    
    % setup roots vector
    roots{k+1}=zeros(k+1,1);
    
    % anonymous function to compute the polynomial... wrapper for
    % fzero
    f = @(x) hval(k+1,x);
    
    for j=1:k+1
        
        % loop through each interval and use fzero to find the root in that
        % interval
        roots{k+1}(j) = fzero(f,[ints(j) ints(j+1)]);
        
    end
    
end

% the quadrature points are the roots of the degree n polynomial
pp=roots{n+1};

% weights come from Wikipedia Gauss Legendre quadrature
% Note: modify for probability measure 0.5 on [-1,1]
hd=(n+1)*( pp.*hval(n+1,pp)-hval(n,pp) )./(pp.^2-1);
pw=1./((1-pp.^2).*hd.^2);

% Note: standard version would be 
% pw=2./((1-pp.^2).*hd.^2);
% but we use probabilists scaling.

%--------------------------------------
% function to compute the Hermite polynomial
% of degree n
%--------------------------------------

function z=hval(n,x)

% this computes the values of polynomials of degrees 0,...,n
tmp=mylegendre(n,x);

% ...we just want the degree n one.
z=tmp(:,end);