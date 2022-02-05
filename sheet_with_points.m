% 2D surface in R^3 with tensor product mesh points
%
%  sp = sheet_with_points(s,ruleu,rulev) represents a surface augmented
%  with tensor product mesh points on the surface. s is a sheet
%  representing the surface and ruleu and rulev are of class rule
%  (representing quadrature rules on [0,1]). The Mesh points are a tensor
%  product of the quadrature points associated with ruleu and rulev.
%
%  sp.points is a 3 x n matrix with sp.points(:,k) storing the x,y,z
%  coordinates of the kth point.
%
%  sp.plot() plots the mesh points.
%
%  h = sp.plot() returns the graphics handle to the mesh points.
%
%  sp.visualise() visualises the surface, and plots the surface normal and
%  the mesh points.
%
%  Note: sheet_with_points is a child class of sheet and has the same methods.
%
% See also: sheet, translated_sheet.
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


classdef sheet_with_points < handle
    
    properties
        sheet
        ruleu
        rulev
        pointsu
        pointsv
        weightsu
        weightsv
        points
        weights
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = sheet_with_points(sheet,ruleu,rulev)
            
            % store the input parameters
            self.sheet = sheet;
            self.ruleu = ruleu;
            self.rulev = rulev;
            
            % generate tensor product of the mesh points
            [self.pointsu,self.pointsv] = meshgrid(ruleu.points(),rulev.points());
            
            % generate tensor product of the corresponding quadrature
            % weights
            [self.weightsu,self.weightsv] = meshgrid(ruleu.points(),rulev.points());

            % combine the ruleu and rulev weights into a single quadrature
            % weight
            self.weights = self.weightsu(:) .* self.weightsv(:);
            
            % calculate the image of the mesh points on the surface and
            % store them
            self.points = sheet.evaluate(self.pointsu(:),self.pointsv(:));
            
        end
        
        %-----------------------------------------
        % plot the mesh points
        %-----------------------------------------

        function varargout = plot(self,opts)
            
            % set default for opts
            if nargin < 2
                opts = 'kx';
            end
            
            % plot the mesh points
            h = plot3(self.points(1,:),self.points(2,:),self.points(3,:),opts);
            
            % if an output is required then return the handle h
            if nargout > 0
                varargout{1} = h;
            end
            
        end
        
        %-----------------------------------------
        % visualise
        %-----------------------------------------

        function visualise(self,color)
            
            if nargin < 2
                color = 'k';
            end
            
            hold_state = ishold;
            
            self.sheet.visualise(color);
            
            hold on
            
            plot3(self.points(1,:),self.points(2,:),self.points(3,:),[color 'x'])
            
            if ~hold_state
                hold off
            end
            
        end
        
    end
    
end