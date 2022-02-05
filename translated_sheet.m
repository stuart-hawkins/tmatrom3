% 2D surface in R^3 parametrised by two parameters
%
%  s1 = translated_sheet(x,s) represents a surface translated by x where x
%  is a 3 x 1 vector. The original surface is represented by the sheet s.
%
%  Note: translated_sheet is a child class of sheet and has the same methods.
%
%  Note: the translation could easily be built into the mapping in the
%  sheet class. However we may want to make use of other child classes of
%  sheet that define particular surfaces eg the surface of a sphere.
%  Translation allows us to apply to translation of such predefined
%  surfaces.
%
% See also: sheet_with_points, sheet.
%
% Stuart C. Hawkins - 13 August 2021

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


classdef translated_sheet < sheet
    
    properties        
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = translated_sheet(translation,sheet_in)

            % define a new parametrisation of the surface based on the one
            % in sheet_in but with the translation added
            f = @(u,v) repmat(translation,1,length(v)) + sheet_in.f(u,v);
            
            % call the parent constructor... note that a constant
            % translation doesn't change the partial derivatives
            self = self@sheet(f,sheet_in.du,sheet_in.dv);
            
        end
        
    end
    
end
        
