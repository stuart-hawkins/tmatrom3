% Dirichlet relation class
%
%  obj = dirichlet_relation(weights,bases,surface,fun) represents a
%  Dirichlet boundary condition applied to bases on a surface represented 
%  by an object of class sheet. In particular,
%
%    weights(1) * bases(1) + ... + weights(end) * bases(end) + fun = 0.
%
%  Here weights is a vector, bases is a vector of class basis, and surface
%  is of class sheet_with_points.
%
%  See relation for more details and an example.
%
% See also: relation, neumann_relation, basis, sheet_with_points.
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


classdef dirichlet_relation < relation
    
    properties
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = dirichlet_relation(weights,bases,sheet,fun)
            
            % call parent constructor
            self = self@relation(weights,bases,sheet,fun);
            
        end
   
        %-----------------------------------------
        % compute scattering matrix
        %-----------------------------------------

        function matrix = inner_matrix(self)

            % matrix involves evaluating self.bases on
            % self.curve
            matrix = self.bases.field(self.sheet,self.weights);
            
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
                
                % then evaluate the incident field
                rhs = self.fun.evaluate(self.sheet.points);
                
            else
                
                % report an error                
                error('fun is not of a supported type')
                
            end

            % make sure the right hand side is a column vector
            rhs = rhs(:);
            
        end
        
    end
    
end