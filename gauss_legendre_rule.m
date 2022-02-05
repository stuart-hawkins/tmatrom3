% Gauss-Legendre quadrature rule class
%
%  q = gauss_legendre_rule(n) creates object q representing a Gauss-Legendre 
%  rule on [0,1] with n points.
%
%  x = q.points() returns vector x containing the quadrature points.
%
%  w = q.weights() returns vector w containing the corresponding quadrature
%  weights.
%
% Example:
%
%    sum(w.*f(x)) approximates the integral of f over the interval [0,1].
%
% Stuart C. Hawkins - 23 August 2021

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


classdef gauss_legendre_rule < rule

    properties
    end
    
    methods
    
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = gauss_legendre_rule(n)
            
            % call parent constructor
            self = self@rule(n);
            
        end
        
        function val = points(self)
           
            % compute the Gauss-Legendre points
            [pp,pw]=gauss_legendre(self.n-1);
            
            % transform from x to theta coordinates... and note we are
            % parametrising from [0,1] rather than [0,pi] so divide by pi
            val = acos(pp(:))/pi;
                        
        end
        
        function val = weights(self)
           
            % compute the Gauss-Legendre points and weigths
            [pp,pw]=gauss_legendre(self.n-1);
            
            % return the weights
            val = pw;
                        
        end
        
    end
    
end