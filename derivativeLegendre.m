% Derivative of the associated Legendre functions
%
%  p = derivativeLegendre(degree,x) instantiates an object to compute
%  derivatives of the associated Legendre functions.
%
%  p.setup(x) computes and stores the values of the derivative of the 
%  associated Legendre function at the points x.
%
%  P = p.get(n) for n<=degree returns the values P of the derivative of the 
%  associated Legendre functions of degree n in a length(x) x (2*n+1) matrix.
%  P(:,m+n+1) holds the values of the polynomial of order m.
%
%  P = p.get(n,m) returns the values of the derivative of the associated 
%  Legendre function of degree n and order m at the points x.
%
% Note: the polynomials are normalised as in Equation (3.2) in
% Ganesh and Hawkins, Numerical Algorithms (2006) 43:25-60.
%
% Note: this code uses derivative_legendre.m downloaded from [1].
%
% See also: associatedLegendre, derivative_legendre.
%
% References: 
%
%  [1] https://www.mathworks.com/matlabcentral/fileexchange/46930-first-derivative-of-normal-associated-legendre-polynomials
%
% Stuart C. Hawkins - 20 April 2021

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


classdef derivativeLegendre < handle
    
    properties
        degree
        data
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = derivativeLegendre(degree,x)
            
            self.degree = degree;
            self.data = [];
            
            self.setup(x);
            
        end
        
        %-----------------------------------------
        % setup
        %-----------------------------------------

        % This computes the derivative of the associated Legendre function 
        % values at the points x and stores them for later use.

        function setup(self,x)
            
            for n = 0:self.degree
                
                % get Schmidt seminormalised polynomial values
                if n>0
                    schmidt = squeeze(legendre_derivative(n,x,'sch')).';
                else
                    schmidt = zeros(2*n+1,length(x));
                end
                
                % setup sparse matrix to implement normalisation as per (3.2) in
                % Ganesh and Hawkins, Numerical Algorithms (2006) 43:25-60
                if n==0
                    lambda = sqrt(1/(4*pi));
                else
                    % negative m scaling
                    scalem = sqrt((2*n+1)/(8*pi)) * ones(n,1);
                    
                    % positive m scaling
                    scalep = sqrt((2*n+1)/(8*pi)) * [sqrt(2);ones(n,1)];
                    
                    % polynomial for m and -m both come from the m column
                    % in Schmidt... lamda is of the form [M,P] where M does
                    % scaling for -m and P does scaling for m.
                    lambda = [fliplr(spdiags(scalem,-1,n+1,n)), ...
                        spdiags((-1).^(0:n).'.*scalep,0,n+1,n+1)];
                end
                
                % post multiply by lambda to implement normalisation
                self.data{n+1} = schmidt * lambda;
                
            end
            
        end
        
        %-----------------------------------------
        % get
        %-----------------------------------------

        function val = get(self,varargin)
            
            if nargin==2
            
                % n = varargin{1}
                val = self.data{varargin{1}+1};
                
            elseif nargin == 3
                
                % n = varargin{1}
                % m = varargin{2}
                
                tmp = self.data{varargin{1}+1};
                
                % Note: the (n,m) entry is in the (n+m+1)th column of
                % tmp
                val = tmp(:,varargin{1}+1+varargin{2});
            
            else
                
                error('get takes two or three parameters')
                
            end
                
        end
        
    end
    
end