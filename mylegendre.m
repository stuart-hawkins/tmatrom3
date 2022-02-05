% Legendre polynomial
%
%    H = mylegendre(n,x) computes the values H of the Legendre polynomials of
%    degrees 0,...,n at points given in x. The values of the degree k
%    polynomial are in H(:,k+1).
%
%    H = mylegendre(n,x,'N') computes the values of the normalized Legendre
%    polynomial.
%
% See also: hermite.
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


function H=mylegendre(n,x,opts)

% set default for opts
if nargin<3
    opts='';
end

% make sure x is a column vector
x=x(:);

% setup return array
m=length(x);
H=zeros(m,n+1);

% set values for degree 0 polynomials
% Note: degree k polynomial is stored in H(:,k+1)
H(:,1) = ones(m,1);

if n>0
    
    % compute degree 1 values
    H(:,2) = x;
    
    for k=1:n-1
        
        % use three term recurrence for higher degree polynomials
        H(:,k+2) = ((2*k+1)/(k+1))*x.*H(:,k+1)-(k/(k+1))*H(:,k);
        
    end
    
end

% normalise the polynomials if necessary
if strfind(opts,'N')
    
    % Note: we assume probability measure on [-1,1] is 0.5 so normalise wrt
    % to this measure
    
    for k=0:n
        
        H(:,k+1)=H(:,k+1)*sqrt((2*k+1));
        
        % Note standard normalised version would be
        % H(:,k+1)=H(:,k+1)*sqrt((2*k+1)/2);
        % but we use probabilists scaling.
        
    end
    
end