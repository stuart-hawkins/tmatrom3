% 2D disk in R^3 parametrised from [0,1] x [0,1]
%
%  s = disk(radius) represents the disk in the xy-plane centred at the
%  origin, with radius r.
%
%  F = s.evaluate(u,v) computes the coordinates of points on the disk.
%  If u and v are n x 1 matrix then F is a 3 x n matrix and F(:,k) holds
%  the (x,y,z) coordinates of the point on the disk corresponding to the
%  point (u(k),v(k)).
%   
%  This is a child class of sheet. See sheet for a full list of methods.
%
% See also: sheet, translated_sheet, spherical_segment, ellipsoid_segment.
%
% Stuart C. Hawkins - 24 August 2021

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


classdef disk < sheet
    
    properties
    end
    
    methods
   
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = disk(radius)
                        
            % set default for radius
            if nargin<1
                radius = 1;
            end
            
            % specify polar coordinates parametrisation of the disk 
            % for the surface mapping....
            f = @(r,theta) radius*[r.*cos(2*pi*theta);r.*sin(2*pi*theta);zeros(size(theta))];
            
            % ....and its derivative with respect to u (r)...
            dr = @(r,theta) radius*[2*pi*cos(theta);2*pi*sin(theta);zeros(size(theta))];
            
            % ....and its derivative with respect to v (theta)...            
            dtheta = @(r,theta) radius*[-2*pi*r.*sin(theta);2*pi*r.*cos(theta);zeros(size(theta))];
            
            % pass the mapping and its partial derivatives to the parent
            % constructor
            self = self@sheet(f,dr,dtheta);
            
        end
        
    end
    
end