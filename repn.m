% Wigner coefficients
%
%  d = repn(dpi2,alpha,n) computes the Wigner coefficients of order n for 
%  rotation of the coordinate system by angle alpha. dpi2 are the
%  coefficients for angle alpha = pi/2. The coefficients for alpha are
%  computed from the coefficients for theta = pi/2 using the formula (3.31)
%  in [1].
%
% See also: repnpi2.
%
% References:
%
% [1] Ganesh and Graham, A high order algorithm for obstacle scattering
% in three dimensions, Journal of Computational Physics 198 (2004)
% 211--242.
%
% Stuart C. Hawkins - 26 October 2021

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


function d = repn(dpi2,alpha,n)

% Note: for n=0 the product with a 1 x 1 spdiags matrix gives a sparse
% matrix as the result. We do n=0 separately to avoid this issue.

if n==0
    
    % get vector of j values
    j = -n:n;
    
    % compute d using the formula (3.31) in [1]
    d = exp((1i*pi/2)*(j-j.')).*(dpi2.' * diag(exp(1i*j.'*alpha)) ...
        * dpi2);

else

    % get vector of j values
    j = -n:n;
    
    % compute d using the formula (3.31) in [1]
    d = exp((1i*pi/2)*(j-j.')).*(dpi2.' * spdiags(exp(1i*j.'*alpha),0,2*n+1,2*n+1) ...
        * dpi2);
    

end