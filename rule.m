% Quadrature rule class
%
%  q = rule(n) creates object q representing a quadrature rule on [0,1] 
%  with n points.
%
%  x = q.points() returns vector x containing the quadrature points.
%
% Note: This is an abstract class and cannot be instantiated.
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


classdef rule < handle

    properties
        n
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = rule(n)
            
            self.n = n;
            
        end
        
    end
    
    %-----------------------------------------
    % abstract function definitions
    %-----------------------------------------
    
    % Note: these specify the interface for methods that must be provided
    % by child classes

    methods(Abstract)
           
        val = points(self);        
        
    end
    
end
        