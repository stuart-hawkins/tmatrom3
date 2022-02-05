% Incident field class
%
%  This abstract class provides the key interfaces for objects of class 
%  incident. These are used by objects of solver class for creating T-matrices.
%
%  Note: this is an abstract class, it cannot be instantiated.
%
% See also: plane_wave, point_source.
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



classdef (Abstract) incident < handle
    
    properties
    end
    
    methods
        
    end
    
    methods(Abstract=true)

        % these must be overloaded in the child class
        cof = get_coefficients(self,centre,nmax);
        val = evaluate(self,points,mask);
                
    end
    
end