% Regular spherical wavefunction class
%
%  e = regularwavefunction(l,j,kwave,x0) returns e, which represents a regular 
%  wavefunction with degree l and order j, wavenumber kwave and origin x0.
%
%  f = e.evaluate(points) evaluates the wavefunction at the points 
%  specified in the 3 x n matrix points.
%
%  [dx,dy,dz] = e.evaluateGradient(points) evaluates the gradient of the 
%  wavefunction (dx,dy,dz) at the points specified in the 3 x n matrix 
%  points.
%
%  G = e.evaluateGradient(points) returns the gradient in a 3 x n matrix.
%
%  cof = e.get_coefficients(x0,n) returns the regular wavefunction
%  expansion coefficients of the wavefunction. The expansion is centred at
%  x0 and has order n.
%
% See also: wavefunction, radiatingwavefunction.
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


classdef regularwavefunction < wavefunction
    
    properties
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = regularwavefunction(varargin)
            
            % call parent constructor
            self = self@wavefunction(varargin{:});
            
        end
       
        %-----------------------------------------
        % the part of the wavefunction that depends
        % on r
        %-----------------------------------------

        % Note: this is used by the parent wavefunction class for
        % evaluating the wavefunction
        
        function val = radial_function(self,r)
           
            val = sphbesselj(self.degree,r);
            
        end
        
        %-----------------------------------------
        % the derivative of the part of the wavefunction 
        % that depends on r
        %-----------------------------------------

        % Note: this is used by the parent wavefunction class for
        % evaluating the gradient of the wavefunction

        function val = derivative_radial_function(self,r)
           
            val = sphbesseljd(self.degree,r);
            
        end
        
        %-----------------------------------------
        % get the expansion coefficients of the wavefunction
        %-----------------------------------------

        function cof = get_coefficients(self,centre,nmax)
    
            % check centre
            if max(abs(centre)) ~= 0                
                error('Only centre=[0;0;0] currently supported.')                
            end

            % setup return vector
            v = zeros((nmax+1)^2,1);

            % convert to a cell
            c = vec2cell(v);
            
            if self.degree <= nmax
                % set entry corresponding to current order to 1
                c{self.degree+1}(self.order+self.degree+1) = 1;
            else
                % the order is larger than nmax and so the coefficient
                % vector does not have a place to put the 1
                warning('order > nmax so coefficient vector is incomplete.')
            end                           
            
            % convert to a vector
            cof = cell2vec(c);

        end
        
    end
    
end