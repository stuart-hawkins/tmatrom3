% Neumann relation class
%
%  obj = neumann_relation(weights,bases,sheet,fun) represents a Neumann
%  boundary condition applied to bases on a surface represented by an object
%  of class sheet. In particular,
%
%    weights(1) * D bases(1) + ... + weights(end) * D bases(end) + fun = 0
%
%  where D denotes the normal derivative to surface.
%
%  See relation for more details and an example.
%
%  Here weights is a vector, bases is a vector of class basis, and surface
%  is of class sheet_with_points.
%
% See also: relation, dirichlet_relation, basis, sheet_with_points.
%
% Stuart C. Hawkins - 21 April 2021

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


classdef neumann_relation < relation
    
    properties
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = neumann_relation(weights,bases,sheet,fun)
            
            % call parent constructor
            self = self@relation(weights,bases,sheet,fun);
            
        end
   
        %-----------------------------------------
        % compute scattering matrix
        %-----------------------------------------

        function matrix = inner_matrix(self)

            % matrix involves evaluating the normal trace of self.bases on
            % self.curve
            matrix = self.bases.normaltrace(self.sheet,self.weights);
                        
        end            
        
        %-----------------------------------------
        % compute right hand side vector
        %-----------------------------------------

        function rhs = inner_rhs(self)
            
            % what we do depends on the type of self.fun....
            if isa(self.fun,'function_handle')
                
                % then evaluate the function
                rhs = self.fun(self.sheet.points);
                
            elseif isa(self.fun,'incident')
                
                % then evaluate the gradient and take its dot product with
                % the normal

                % get the normal
                n = self.sheet.sheet.normal(self.sheet.pointsu,self.sheet.pointsv);
                
                % get the gradient of the incident wave
                [dx,dy,dz] = self.fun.evaluateGradient(self.sheet.points);
                
                % compute the dot product
                rhs = n(1,:).*dx + n(2,:).*dy + n(3,:).*dz;
            
            else
                
                % report an error
                error('fun is not of a supported type')
                
            end
            
            % make sure the right hand side is a column vector
            rhs = rhs(:);
            
        end
        
    end
    
end