% Segment of a sphere in R^3 parametrised from [0,1] x [0,1]
%
%  s = spherical_segment(Itheta,Iphi,radius) represents a segment of the
%  sphere of radius r centred at the origin parametrised in spherical polar
%  coordinates
%
%    x = r * sin(theta) * cos(phi)
%    y = r * sin(theta) * sin(phi)
%    z = r * cos(theta)
%
%  The segment is that part of the sphere with theta between Itheta(1) and
%  Itheta(2), and phi between Iphi(1) and Iphi(2). If Itheta = [0,pi] and
%  Iphi = [0,2*pi] then the segment is the whole sphere.
%
%  F = s.evaluate(u,v) computes the coordinates of points on the segment.
%  If u and v are n x 1 matrix then F is a 3 x n matrix and F(:,k) holds
%  the (x,y,z) coordinates of the point on the segment corresponding to the
%  point (u(k),v(k)).
%   
%  This is a child class of sheet. See sheet for a full list of methods.
%
% See also: sheet, translated_sheet, ellipsoid_segment, disk.
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


classdef spherical_segment < sheet
    
    properties
    end
    
    methods
   
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = spherical_segment(thetalim,philim,radius)
            
            % set default for radius
            if nargin<3
                radius = 1;
            end
            
            % map u in [0,1] to theta in thetalim
            theta = @(u) thetalim(1) + (thetalim(2)-thetalim(1)) * u;

            % map v in [0,1] to phi in philim
            phi = @(v) philim(1) + (philim(2)-philim(1)) * v;

            % ... and the derivatives of those mappings wrt u and v
            thetad = (thetalim(2)-thetalim(1));
            phid = (philim(2)-philim(1));

            % specify spherical polar coordinates parametrisation of the sphere 
            % for the surface mapping....
            f = @(u,v) radius*[sin(theta(u)).*cos(phi(v));sin(theta(u)).*sin(phi(v));cos(theta(u))];
            
            % ....and its derivative with respect to u (theta)...
            dtheta = @(u,v) radius*[thetad*cos(theta(u)).*cos(phi(v));thetad*cos(theta(u)).*sin(phi(v));-thetad*sin(theta(u))];

            % ....and its derivative with respect to v (phi)...
            dphi = @(u,v) radius*[-phid*sin(theta(u)).*sin(phi(v));phid*sin(theta(u)).*cos(phi(v));zeros(size(theta(u)))];
            
            % pass the mapping and its partial derivatives to the parent
            % constructor
            self = self@sheet(f,dtheta,dphi);
            
        end
        
    end
    
end