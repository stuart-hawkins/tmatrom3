% Plane wave incident field.
%
%  u = plane_wave(d,k) returns a plane wave object u with 
%  wavenumber k and direction given by the unit vector d (with length 3).
%
% Also:
%
%   f = u.evaluate(p) returns the values f of the plane wave at points p.
%   Here p must be a 3 x n matrix.
%
%   f = u.evaluate(z,mask) returns the values f of the plane wave at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy,dz] = u.evaluateGradient(p) returns dx, dy and dz the partial 
%   derivatives of the plane wave in the x, y and z directions respectively
%   at the points p. Here p must be a 3 x n matrix.
%
%   [dx,dy,dz] = u.evaluateGradient(z,mask) returns dx, dy and dz the partial 
%   derivatives of the plane wave in the x, y and z directions respectively
%   at the points p for which mask==1 and NaN elsewhere.
%
%   cof = u.get_coefficients(x0,n) returns the vector cof of regular
%   wavefunction expansion coefficients of the plane wave field with 
%   wavefunction origin x0 and order n.
%
% See also: point_source, incident.
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


classdef plane_wave < incident
    
    properties
        direction
        kwave
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = plane_wave(direction,kwave)
            
            % check that direction is a unit vector
            if abs(norm(direction)-1) > 1e-15
                direction = direction/norm(direction);
                warning('direction was not a unit vector and has been rescaled')
            end
            
            % set wavenumber
            self.kwave = kwave;
            
            % set direction (as column vector)
            self.direction = direction(:);
            
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
            
            % evaluate incident field
            v = exp(1i*self.kwave*self.direction(:).' * points);
            
            % insert values into the return array
            if nargin>2
                val(:,mask) = v;
            else
                val = v;
            end
            
            % reshape the return array to match points
            if length(n)==2
                val = reshape(val,n(end),1);
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
            
            % evaluate incident field
            v = exp(1i*self.kwave*self.direction(:).' * points);
            
            % insert values into the return array
            if nargin>2
                dx(:,mask) = 1i*self.kwave*self.direction(1)*v;
                dy(:,mask) = 1i*self.kwave*self.direction(2)*v;
                dz(:,mask) = 1i*self.kwave*self.direction(3)*v;
            else
                dx = 1i*self.kwave*self.direction(1)*v;
                dy = 1i*self.kwave*self.direction(2)*v;
                dz = 1i*self.kwave*self.direction(3)*v;
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
        % get coefficients
        %-----------------------------------------

        function cof = get_coefficients(self,centre,nmax)
        
            % get polar coordinates for the wave direction
            theta = acos(self.direction(3));
            phi = atan2(self.direction(2),self.direction(1));
            
            % evaluate the spherical harmonic for the wave direction
            Y = associatedLegendre(nmax,cos(theta));
            
            % set up cell array holding indexes... cof{n+1} holds the
            % coefficients for order n
            for n=0:nmax

                j = -n:n;
                
                % set the coefficients
                cof{n+1} = 4*pi*1i^n * Y.get(n) .* exp(-1i*j*phi);
                
            end
            
            % convert the cell array into a vector
            cof = cell2vec(cof);

            % adjust the coefficients if the centre of the scatterer is not
            % the origin
            if max(abs(centre))~=0
                dp = dot(self.direction,centre);
                phase = exp(1i*self.kwave*dp);
                cof = phase * cof;
            end

        end
        
    end
    
end
       