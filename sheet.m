% 2D surface in R^3 parametrised by two parameters
%
%  s = sheet(f,du,dv) represents a surface parametrised by a function
%  f(u,v) for u,v in [0,1]. The function f must map n x 1 matrices u and v 
%  into a 3 x n matrix. The functions du and dv specify the partial
%  derivatives of f with respect to the parameters u and v respectively.
%  The functions du(u,v) and dv(u,v) map n x 1 matrices into 3 x n matrices.
%
%  F = s.evaluate(u,v) computes the coordinates of points on the surface.
%  If u and v are n x 1 matrix then F is a 3 x n matrix and F(:,k) holds
%  the (x,y,z) coordinates of the point on the surface corresponding to the
%  point (u(k),v(k)).
%
%  [Du,Dv] = s.derivative(u,v) computes the partial derivatives of the
%  mapping representing the surface with respect to the parameters u and v.
%  If u and v are n x 1 matrices then Du (resp. Dv) is a 3 x n matrix and Du(:,k) 
%  (resp. Dv(:,k)) holds the (x,y,z) coordinates of the partial derivative of the
%  surface parametrisation with respect to u (resp. v) at the point on the surface
%  corresponding to the point (u(k),v(k)).
%  
%  J = s.jacobian(u,v) computes the jacobian of the mapping representing 
%  the surface with respect to the parameters u and v. If u and v are n x 1 
%  matrices then J is an n x 1 matrix and J(k) holds the (x,y,z) coordinates of 
%  the Jacobian of the surface parametrisation at the point on the surface
%  corresponding to the point (u(k),v(k)).
%
%  N = s.normal(u,v) computes the normal to surface. If u and v are n x 1 
%  matrices then N is a 3 x n matrix and N(:,k) holds the (x,y,z) coordinates 
%  of the normal to the surface at the point on the surface corresponding 
%  to the point (u(k),v(k)).
%
%  s.plot() visualises the surface using a surface plot of the surface and
%  the Jacobian.
%
%  h = s.plot() returns the graphics handle to the surface plot.
%
%  s.visualise() visualises the surface and plots the surface normal.
%
% See also: sheet_with_points, translated_sheet.
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


classdef sheet < handle
    
    properties
       f
       du
       dv
    end
    
    methods
   
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = sheet(f,du,dv)
            
            % store functions
            self.f = f;
            self.du = du;
            self.dv = dv;
            
        end
        
        %-----------------------------------------
        % compute coordinates of the surface
        %-----------------------------------------

        function val = evaluate(self,u,v)
           
            val = self.f(u(:).',v(:).');
            
        end
        
        %-----------------------------------------
        % compute the partial derivatives of the surface
        % with respect to the parameters u and v
        %-----------------------------------------

        function [du,dv] = derivative(self,u,v)
            
            du = self.du(u(:).',v(:).');
            dv = self.dv(u(:).',v(:).');
            
        end
                
        %-----------------------------------------
        % compute the jacobian of the transformation
        % from (u,v) coordinates to the surface
        % in (x,y,z) coordinates
        %-----------------------------------------

        function val = jacobian(self,u,v)
            
            % get the partial derivative wrt u and v
            [du,dv] = self.derivative(u(:).',v(:).');
            
            % compute the jacobian as the norm of the cross product of the
            % partial derivatives
            val = sqrt(sum(cross(du,dv,1).^2,1));
            
        end
        
        %-----------------------------------------
        % compute the  normal to the surface
        %-----------------------------------------

        function val = normal(self,u,v)
            
            % get the partial derivatives wrt u and v
            [du,dv] = self.derivative(u(:).',v(:).');
            
            % the normal points in the direction of the cross product...
            val = cross(du,dv,1);
            
            % compute the norm of the cross product
            nrm = sqrt(sum(val.^2,1));
            
            % divide by the norm to get the unit normal
            val = val ./ repmat(nrm,3,1);
            
        end

        %-----------------------------------------
        % visualise the surface and the jacobian
        %-----------------------------------------

        function varargout = plot(self)

            % number of mesh points used to represent the surface
            n = 30;
            
            % set up meshes in the u and v coordinates
            u = linspace(0,1,n);
            v = linspace(0,1,n);
            
            % turn the u and v meshes into 2D grid points
            [uu,vv] = meshgrid(u,v);
            
            % evaluate the surface coordinates at the grid points
            f = self.evaluate(uu(:),vv(:));
            
            % now f is a 3 x n^2 matrix... turn it into a 3 x n x n array
            f = reshape(f,3,n,n);
            
            % get the jacobian at the mesh points
            J = self.jacobian(uu(:),vv(:));
            
            % now J is a n^2 vector... turn it into an n x n array
            J = reshape(J,n,n);
            
            % do a surf plot of the surface
            h = surf(squeeze(f(1,:,:)),squeeze(f(2,:,:)),squeeze(f(3,:,:)),J);
            
            % if an output is required then return the handle for the surf
            % plot
            if nargout > 0
                varargout{1} = h;
            end
            
        end
        
        %-----------------------------------------
        % visualise the surface
        %-----------------------------------------

        function visualise(self)
                        
            % record the hold state so that we can restore it later
            hold_state = ishold;
            
            % number of mesh points used to represent the surface
            n = 60;
            
            % set up meshes in the u and v coordinates
            u = linspace(0,1,n);
            v = linspace(0,1,n);
            
            % turn the u and v meshes into 2D grid points
            [uu,vv] = meshgrid(u,v);
            
            % evaluate the surface coordinates at the grid points
            f = self.evaluate(uu(:),vv(:));
            
            % now f is a 3 x n^2 matrix... turn it into a 3 x n x n array
            f = reshape(f,3,n,n);
            
            % get the jacobian at the mesh points
            J = self.jacobian(uu(:),vv(:));
            
            % now J is a n^2 vector... turn it into an n x n array
            J = reshape(J,n,n);
            
            % do a surf plot of the surface
            surf(squeeze(f(1,:,:)),squeeze(f(2,:,:)),squeeze(f(3,:,:)),J);
                 
            % hold the figure ready to add the normals
            hold on
            
            % number of mesh points used to visualise the normals
            m = 15;
            
            % set up meshes in the u and v coordinates
            u1 = linspace(0,1,m);
            v1 = linspace(0,1,m);
            
            % turn the u and v meshes into 2D grid points
            [uu1,vv1] = meshgrid(u1,v1);
            
            % evaluate the surface coordinates at the grid points
            f1 = self.evaluate(uu1(:),vv1(:));
            
            % get the normal at the mesh points
            nrml = self.normal(uu1(:),vv1(:));
            
            % do a quiver plot of the normals
            quiver3(f1(1,:),f1(2,:),f1(3,:),nrml(1,:),nrml(2,:),nrml(3,:))
            
            % restore the hold state
            if ~hold_state
                hold off
            end            
            
        end
        
    end
    
end