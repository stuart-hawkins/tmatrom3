% Relation class
%
%  obj = relation(weights,bases,sheet,fun) creates an object of class
%  relation. 
%
%  A = obj.matrix() returns a matrix A with block structure 
%  A = [A_1,A_2,...,A_N] when bases is an 1 x N array. The entry A_k(i,j) 
%  gives the field from the jth source in bases(k) at the ith point in
%  sheet. The "kind" of field depends on the type of relation.
%
%  When obj is an M x 1 relation array the matrix A has block structure
%  A = [A_11,...,A_1N;A_21,...,A_2N;...;A_M1,...,A_MN]. The entry A_lk(i,j)
%  gives the field from the jth source in bases(k) at the ith point. The 
%  "kind" of field depends on the type of the relation obj(l).
%
%  f = obj.rhs() returns a matrix f with block structure f=[f_1;...;f_M]
%  when obj is an M x 1 relation array. The entry f_l(i) gives the incident
%  field at the ith point of sheet. The "kind" of field depends on the type
%  of the relation obj(l).
%
%  This class exists just to handle vectors of relations, which are used to
%  specify multiple boundary conditions between bases and sheets. The
%  objects in the vector will be instances of the child classes of
%  relation, eg dirichlet_relation and neumann_relation.
%
%  Example:
%
%    Apply a transmission boundary condition across a boundary represented
%    by a sheet S. Letting u and v be bases generating exterior and interior 
%    fields respectively, and w an incident field, the boundary conditions 
%    are
%
%       -u + v = w
%       -Du + Dv = Dw
%
%    where D denotes the normal derivative.
%
%    This boundary condition is represented as
%
%     rel = 
%      [
%         dirichlet_relation([-1 1],[u v],S,w);
%         neumann_relation([-1 1],[u v],S,w)]
%      ]
%
%    Now
%
%     cof = rel.matrix() \ rel.rhs()
%
%    gives the weights for bases u and v that satisfy the boundary
%    conditions in the vector cof.
%
%  Note: this is an abstract class, it cannot be instantiated.
%
% See also: solver, dirichlet_relation, neumann_relation.
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


% Note: we subclass matlab.mixin.Heterogeneous so that we can handle arrays
% of relations.

classdef relation < handle & matlab.mixin.Heterogeneous
    
    properties
        weights
        bases
        fun
        sheet
    end
    
    methods
   
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = relation(weights,bases,sheet,fun)
        
            % store parameters
            self.weights = weights;
            self.bases = bases;
            self.fun = fun;
            self.sheet = sheet;
            
        end
        
        
    end
    
    methods(Sealed)
        
        %-----------------------------------------
        % construct scattering matrix
        %-----------------------------------------

        function matrix = matrix(self)
           
            % initialise output
            matrix = [];

            % loop through vector of relations
            for j=1:length(self)
                
                % self(j) will be an instance of one of the child classes
                % of relation... call that child class's inner_matrix()
                % method to assemble that part of the matrix
                matrix = [matrix;self(j).inner_matrix()];
                
            end
            
        end
        
        %-----------------------------------------
        % assemble right hand side
        %-----------------------------------------

        function rhs = rhs(self)
           
            % initialise output
            rhs = [];
            
            % loop through vector of relations
            for j=1:length(self)
                
                % self(j) will be an instance of one of the child classes
                % of relation... call that child class's inner_rhs()
                % method to assemble that part of the right hand side
                rhs = [rhs;self(j).inner_rhs()];
                
            end
            
        end
        
    end
    
    %-----------------------------------------
    % abstract function definitions
    %-----------------------------------------
    
    % Note: these specify the interface for methods that must be provided
    % by child classes
    
    methods(Abstract)
        
        val = inner_matrix(self);
        val = inner_rhs(self);
        
    end
    
end
    
    