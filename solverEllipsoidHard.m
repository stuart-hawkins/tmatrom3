% Solve acoustic scattering for a sound hard ellipsoid
%
%  s = solverEllipsoidHard(k,inc,x0,ar,n) initialises a solver class to
%  simulate acoustic scattering by an ellipsoid using the method of
%  fundamental solutions. Parameters are:
%
%    k - the wavenumber
%    inc - specifies the incident wave (of class incident or empty)
%    x0 - the centre of the ellipsoid (a 3 x 1 vector)
%    ar - the radius in the x,y,z directions (a 3 x 1 vector)
%    n - parameter controlling the number of discrete sources in the interior
%
%  s = solverEllipsoidHard(k,inc,x0,ar,n,m) further specifies:
%
%    m - parameter controlling the number of surface points at which the
%    boundary condition is matched.
%
% Also:
%
%  s.setup() prepares the class. This only needs to be done once and is
%  independent of the incident wave.
%
%  s.solve() solves the scattering problem. The object s can be
%  subsequently used to compute the far field, cross section etc.
%
%  s.visualise() visualises the scatterer.
%
%  s.setIncidentField(inc) sets the incident field as specified in the cell
%  array inc.
%
%  val = s.getFarField(pts) computes the far field for the
%  incident field inc{1} at points on the unit sphere specified in a 3 x n 
%  matrix pts.
%
%  val = s.getFarField(pts,index) computes the far field for the
%  incident fields inc{index} at points on the unit sphere specified in a 
%  3 x n matrix pts.
%
%  val = s.getCrossSection(pts) computes the cross section for the
%  incident field inc{1} at points on the unit sphere specified in a 3 x n 
%  matrix pts.
%
%  val = s.getCrossSection(pts,index) computes the cross section for the
%  incident fields inc{index} at points on the unit sphere specified in a 
%  3 x n matrix pts.
%  
% See also: solver, solverEllipsoidSoft, solverEllipsoidPenetrable.
%
% Stuart C. Hawkins - 21 April 2021

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


classdef solverEllipsoidHard < solver

    properties
        origin
        shell
        basis_minus
        matrix
        boundary_conditions
        cof
        basis_param
        ar
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solverEllipsoidHard(kwave,incidentField,origin,ar,basis_param,m)

            % set default for m
            if nargin<5
                m = basis_param;
            end
            
            % call parent constructor
            self = self@solver(kwave,incidentField);
            
            % store aspect ratio
            self.ar = ar;

            % store origin
            self.origin = origin;

            % store numnber of sources parameter
            self.basis_param = basis_param;
            
            % setup boundary
            bdry = translated_sheet(self.origin(:),ellipsoid_segment([0 pi],[0 2*pi],ar(1),ar(2),ar(3)));
            self.shell = sheet_with_points(bdry,midpoint_rule(m),midpoint_rule(2*m));
            
            % basis for exterior
            int = translated_sheet(self.origin(:),ellipsoid_segment([0 pi],[0 2*pi],0.95*ar(1),0.95*ar(2),0.95*ar(3)));
            minus = sheet_with_points(int,midpoint_rule(basis_param),rectangle_left_rule(2*basis_param));
            self.basis_minus = basis(minus,kwave);
            
            % record boundary conditions
            self.boundary_conditions = ...
                [neumann_relation(-1,self.basis_minus,self.shell,[])];
            
        end
                 
        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        function setup(self)                      
            
            % compute the system matrix
            self.matrix = self.boundary_conditions.matrix();
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        function solve(self)

            % loop through the right hand sides to assemble the right hand
            % side matrix
            for j=1:length(self.incidentField)
            
                % associate jth incident field with the boundary condition
                self.boundary_conditions(1).fun = self.incidentField{j};
            
                % evaluate the vector associated with the incident field
                rhs(:,j) = self.boundary_conditions.rhs();
                
            end
            
            % Matlab sometimes warns about a rank deficient matrix but this
            % can be safely ignored so turn it off
            warning('off','MATLAB:rankDeficientMatrix')

            % solve the linear system
            % Note: self.matrix is assembled in the setup() method
            self.cof = self.matrix \ rhs;
            
            % turn the warning back on
            warning('on','MATLAB:rankDeficientMatrix')

        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        function val = getFarField(self,points,index)

            % set default for index
            if nargin<3
                index = 1:length(self.incidentField);
            end
            
            % loop through index
            for j=1:length(index)
                
                % compute the far field associated with right hand side
                % index(j)
                val(:,j) = self.basis_minus.farfield(points,self.cof(:,j));
                
            end            

        end
        
        %-----------------------------------------
        % text description of the solver
        %-----------------------------------------
        
        function str = description(self)
           
            str = sprintf('Method of Fundamental Solution solver for sound hard ellipsoid with semi-axes (%0.1f,%0.1f,%0.1f).',...
                self.ar(1),self.ar(2),self.ar(3));
            
        end
        
        %===============================================================
        % you may provide other methods required to implement your solver
        % or help the user
        %===============================================================

        %-----------------------------------------
        % visualise the scatterer
        %-----------------------------------------

        function varargout = visualise(self)
        
            % number of points in each dimension in the surface mesh
            n = 60;
            
            % setup mesh points in the reference rectangle
            u = linspace(0,1,n);
            v = linspace(0,1,n);
            [uu,vv] = meshgrid(u,v);

            % transform the mesh points to the surface
            % Note: sheet.evaluate() method wants points in vectors so we
            % turn uu and vv into column vectors first
            f = self.shell.sheet.evaluate(uu(:),vv(:));
            
            % reshape the result into an n x n matric, corresponding to the
            % grid
            f = reshape(f,3,n,n);
                        
            % surf the surface
            shell = surf(squeeze(f(1,:,:)),squeeze(f(2,:,:)),squeeze(f(3,:,:)));
            
            % tweak the appearance
            set(shell,'facecolor',0.5*[1 1 1])
            set(shell,'edgecolor','none')            

            % return handle of surface if required
            if nargout > 0
                varargout{1} = shell;
            end
            
        end
        
    end % end methods
    
end