% Wigner coefficients
%
%  obj = wigner(order) sets up an object to compute the Wigner
%  coefficients. The coefficients for polar angle rotation theta are
%  computed from the coefficients for theta = pi/2 using the formula (3.31)
%  in [1].
%
%  W = obj.get(theta) gets the coefficients for angle theta. The coefficent
%  d^l_{j,jd} is stored in W{l+1}(j,jd).
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


classdef wigner < handle
    
    properties
        order
        dpi2
    end
    
    methods
   
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = wigner(order)
            
            self.order = order;
           
            self.setup();
            
        end
        
        %-----------------------------------------
        % setup the d(pi/2) matrices
        %-----------------------------------------

        function setup(self)
            
            % setup a cell
            self.dpi2 = cell(self.order+1,1);

            for n=0:self.order
                
                % compute d(pi/2)
                self.dpi2{n+1} = repnpi2(n);
                
            end
            
        end
        
        %-----------------------------------------
        % comopute d(alpha)
        %-----------------------------------------

        function val = get(self,alpha)
            
            % setup a cell
            val = cell(self.order+1,1);
            
            for n=0:self.order
                
                % compute d from d(pi/2)
                val{n+1} = repn(self.dpi2{n+1},alpha,n);
                
            end
            
        end
        
    end
    
end

