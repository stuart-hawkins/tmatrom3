% Point source incident field.
%
%  u = point_source(x0,k) returns a point source object u with 
%  wavenumber k and source location x0 (a vector with length 3).
%
% Also:
%
%   f = u.evaluate(p) returns the values f of the point source at points p.
%   Here p must be a 3 x n matrix.
%
%   f = u.evaluate(z,mask) returns the values f of the point source at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy,dz] = u.evaluateGradient(p) returns dx, dy and dz the partial 
%   derivatives of the point source in the x, y and z directions respectively
%   at the points p. Here p must be a 3 x n matrix.
%
%   [dx,dy,dz] = u.evaluateGradient(z,mask) returns dx, dy and dz the partial 
%   derivatives of the point source in the x, y and z directions respectively
%   at the points p for which mask==1 and NaN elsewhere.
%
%   cof = u.get_coefficients(x0,n) returns the vector cof of regular
%   wavefunction expansion coefficients of the point source field with 
%   wavefunction origin x0 and order n.
%
% See also: plane_wave, incident.
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


classdef point_source < incident
    
    properties
        location
        kwave
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = point_source(location,kwave)           
            
            % set wavenumber
            self.kwave = kwave;
            
            % set location (as column vector)
            self.location = location(:);
            
        end
        
        %-----------------------------------------
        % evaluate
        %-----------------------------------------

        function val = evaluate(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
            
            % initialise return array
            val = zeros(1,size(points,2));
            
            % apply mask if necessary
            if nargin > 2
                points = points(:,mask);
            end
            
            % get distance of points from source
            r = sqrt(sum((points - repmat(self.location,1,size(points,2))).^2,1));
            
            % evaluate incident field
            v = (0.25/pi)*exp(1i*self.kwave*r)./r;
            
            % insert values into the return array
            if nargin>2
                val(:,mask) = v;
            else
                val = v;
            end
            
            % reshape the return array to match points
            if length(n)==2
                val = reshape(val,1,n(end));
            else
                val = reshape(val,n(2:end));
            end
            
        end
        
        %-----------------------------------------
        % evaluate gradient
        %-----------------------------------------

        function [dx,dy,dz] = evaluateGradient(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
            
            % initialise return array
            dx = zeros(1,size(points,2));
            dy = zeros(1,size(points,2));
            dz = zeros(1,size(points,2));
            
            % apply mask if necessary
            if nargin > 2
                points = points(:,mask);
            end
            
            % get direction of points from source
            d = points - repmat(self.location,1,size(points,2));
            
            % get distance of points from source
            r = sqrt(sum(d.^2,1));

            % evaluate incident field
            v = (0.25/pi)*exp(1i*self.kwave*r)...
                .*((1i*self.kwave)*r-1)./r.^3;
            
            % insert values into the return array
            if nargin>2
                dx(:,mask) = v.*d(1,:);
                dy(:,mask) = v.*d(2,:);
                dz(:,mask) = v.*d(3,:);
            else
                dx = v.*d(1,:);
                dy = v.*d(2,:);
                dz = v.*d(3,:);
            end
            
            % reshape the return array to match points
            if length(n)==2
                dx = reshape(dx,1,n(end));
                dy = reshape(dy,1,n(end));
                dz = reshape(dz,1,n(end));
            else
                dx = reshape(dx,n(2:end));
                dy = reshape(dy,n(2:end));
                dz = reshape(dz,n(2:end));
            end
            
        end
        
        %-----------------------------------------
        % evaluate far field
        %-----------------------------------------

        function val = evaluateFarField(self,points)
        
            % compute the far field using (2.15) in Colton and Kress,
            % Inverse Acoustic and Electromagnetic Scattering Theory, 4th
            % edition.
            dp = sum(repmat(self.location,1,size(points,2)) .* points,1);            
            val = (0.25/pi)*exp(-1i*self.kwave*dp);
            
        end
        
        
        %-----------------------------------------
        % get coefficients
        %-----------------------------------------

        function cof = get_coefficients(self,centre,nmax)

            error('not implemented yet')
            
        end
        
    end
    
end
       