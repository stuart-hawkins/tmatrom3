% Midpoint quadrature rule class
%
%  q = midpoint_rule(n) creates object q representing a midpoint quadrature 
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
% See also: rectangle_rule_right, rectangle_rule_left.
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


classdef midpoint_rule < rule

    properties
    end
    
    methods
    
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = midpoint_rule(n)
            
            % call parent constructor
            self = self@rule(n);
            
        end
        
        %-----------------------------------------
        % return points
        %-----------------------------------------

        function val = points(self)
           
            % split [0,1] into n intervals and the points are the midpoints
            % of those intervals
            val = (0.5+(0:self.n-1).')/self.n;
            
        end
        
        %-----------------------------------------
        % return quadrature weights
        %-----------------------------------------

        function val = weights(self)
           
            % return the weights
            val = ones(n,1)/n;
            
        end
        
    end
    
end