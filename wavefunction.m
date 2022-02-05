% Spherical wavefunction class
%
%  e = wavefunction(l,j,kwave,x0) returns e, which represents a wavefunction 
%  with degree l and order j, wavenumber kwave and origin x0.
%
%  f = e.evaluate(points) evaluates the wavefunction at the points 
%  specified in the 3 x n matrix points.
%
%  [dx,dy,dz] = e.evaluateGradient(points) evaluates the gradient of the 
%  wavefunction (dx,dy,dz) at the points specified in the 3 x n matrix 
%  points.
%
%  G = e.evaluateGradient(points) returns the gradient in a 3 x n matrix.
%
%  Note: this is an abstract class, it cannot be instantiated. It
%  implements common methods for the child classes regularwavefunction and
%  radiating wavefunction.
%
% See also: regularwavefunction, radiatingwavefunction.
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


classdef (Abstract) wavefunction < incident & tmatrom3
    
    properties
        degree
        order
        kwave
        origin
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = wavefunction(degree,order,kwave,origin)
            
            % set default for origin
            if nargin < 4
                origin = [0;0;0];
            end
            
            % set properties
            self.degree = degree;
            self.order = order;
            self.kwave = kwave;
            self.origin = origin;
            
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
            
            % subtract origin
            points = points - repmat(self.origin(:),1,size(points,2));
            
            % get radius
            r = sqrt(sum(points.^2,1));
            
            % get azimuth
            phi = atan2(points(2,:),points(1,:));
            
            % evaluate spherical harmonic
            Y = associatedLegendre(self.degree,points(3,:)./r);
            
            % evaluate incident field
            v = self.radial_function(self.kwave * r(:)) .* Y.get(self.degree,self.order) ...
                .* exp(1i*self.order*phi(:));
            
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
        % compute gradient
        %-----------------------------------------

        function varargout = evaluateGradient(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
            
            % initialise return arrays
            dx = zeros(1,size(points,2));
            dy = zeros(1,size(points,2));
            dz = zeros(1,size(points,2));
            
            % apply mask if necessary
            if nargin > 2
                points = points(:,mask);
            end
            
            % subtract origin
            points = points - repmat(self.origin(:),1,size(points,2));
            
            % get radius
            r = sqrt(sum(points.^2,1));
            
            % get azimuth
            phi = atan2(points(2,:),points(1,:));
            theta = acos(points(3,:)./r);
            
            % evaluate spherical harmonic
            Y = associatedLegendre(self.degree,points(3,:)./r);
            YD = derivativeLegendre(self.degree,points(3,:)./r);
            
            % evaluate incident field
            rr = self.radial_function(self.kwave * r(:));
            pp = exp(1i*self.order*phi(:));
            tt = Y.get(self.degree,self.order);
            
            % first we compute the partial derivatives in spherical polar
            % coordinates, then transform to Cartesian coordinates

            % compute partial derivative with respect to r
            vr = self.kwave * self.derivative_radial_function(self.kwave * r(:)) ...
                .* tt .* pp;

            % compute partial derivative with respect to theta
            vt = -sin(theta(:)).*YD.get(self.degree,self.order).* rr .* pp;

            % compute partial derivative with respect to phi
            vp = 1i*self.order * rr .* tt .* pp;
                        
            % compute partial derivative with respect to x
            dx = sin(theta(:)).*cos(phi(:)).*vr ...
                + cos(theta(:)).*cos(phi(:)).*vt./r(:) ...
                - sin(phi(:)).*vp./(sin(theta(:)).*r(:));
            
            % compute partial derivative with respect to y
            dy = sin(theta(:)).*sin(phi(:)).*vr ...
                + cos(theta(:)).*sin(phi(:)).*vt./r(:) ...
                + cos(phi(:)).*vp./(sin(theta(:)).*r(:));

            % compute partial derivative with respect to z
            dz = cos(theta(:)).*vr ...
                - sin(theta(:)).*vt./r(:); 

            % insert values into the return array
            if nargin>2
                dx(:,mask) = dx;
                dy(:,mask) = dy;
                dz(:,mask) = dz;
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
            
            % can return result as a single array or as three components
            if nargout==1
                varargout{1} = [dx;dy;dz];
            elseif nargout==3
                varargout{1} = dx;
                varargout{2} = dy;
                varargout{3} = dz;
            else
                error('Unexpected nargout')
            end
            
        end
        
    end
    
    methods(Abstract=true)
        
        % this must be defined in the child class
        val = radial_function(self,r);
        
    end
    
end