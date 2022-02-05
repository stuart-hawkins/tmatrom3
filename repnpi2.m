% Wigner coefficients
%
%  d = repn(n) computes the Wigner coefficients of order n for 
%  rotation of the coordinate system by angle pi/2. The coefficients are
%  given by (3.26) in [1].
%
% See also: repn.
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


function d = repnpi2(n)

% set the rotation angle to pi/2
alpha = pi/2;

% d_{kd,k} has structure like this:
%
%   ----
%   |\/| A B C
%   ----
%   |/\| D E F
%   ----

% compute region E using (3.32) in [1]
for kd = 0:n    
    for k = -kd:kd
        
        d(kd+n+1,k+n+1) = rep_val(n,abs(kd),abs(k)) ...
            * (0.5)^kd * jacobi0(n-kd,kd-k,kd+k);
    
    end
end

% ...use symmetry relations (3.27) for the other regions
if n>0

    % region F
    for kd = 0:n        
        for k = kd+1:n
            
            d(kd+n+1,k+n+1) = (-1)^(kd-k) * d(k+n+1,kd+n+1);
            
        end        
    end
    
    % regions A and B
    for kd = -n:-1        
        for k = -n:-kd
            
            d(kd+n+1,k+n+1) = (-1)^(kd-k) * d(-kd+n+1,-k+n+1);
            
        end        
    end
    
    % region C
    for kd = -n:-1        
        for k = -kd+1:n
            
            d(kd+n+1,k+n+1) = (-1)^(kd-k) * d(k+n+1,kd+n+1);
            
        end        
    end
    
    % region A
    for kd = 0:n        
        for k = -n:-kd+1
            
            d(kd+n+1,k+n+1) = (-1)^(kd-k) * d(-kd+n+1,-k+n+1);
            
        end        
    end

end

%-----------------------------------------
% compute Jacobi polynomial at 0
%-----------------------------------------

function val = jacobi0(n,a,b)

% setup vector of values
P = zeros(n+1,1);

% setup value of first two polynomials
P(1) = 1;
P(2) = 0.5*(a-b);

% compute remaining polynomials using three term recurrence
for l = 2:n
    P(l+1) = (  (2*l+a+b-1)*( a^2-b^2 ) * P(l) ...
        - 2*(l+a-1)*(l+b-1)*(2*l+a+b) * P(l-1)  ) ...
        /( 2*l*(l+a+b)*(2*l+a+b-2) );
end
   
% return the final value (the highest order polynomial)
val = P(n+1);

%-----------------------------------------
% compute factorial coefficient in (3.32) in [1]
%-----------------------------------------

function val = rep_val(n,j,k)

if j==n && k==n
    val = 1;
elseif k==n
    val = sqrt(2*n+1)*spec_fact(k+1,j-1,n);
else
    val = spec_fact(k,j-1,n);
end

%-----------------------------------------
% compute factorial term in (3.32) in [1]
%-----------------------------------------

function val = spec_fact(j,k,n)

% get integers from j to k
jj = j:k;

% compute double factorial term in one go 
val = prod(sqrt((n+1+jj)./(n-jj)));