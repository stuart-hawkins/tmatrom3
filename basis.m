% Basis of discrete sources on a sheet_with_points.
%
%  B = basis(s,kwave) represents discrete sources at the points of s and
%  wavenumber kwave. Here s must be a sheet_with_points.
%
%  F = field(B,d,w) computes the field from discrete sources associated
%  with bases B multiplied by weights w. Here B is an n vector of class
%  basis and w is an n vector of scalar weights. F(i,j) is the field at
%  point i from source j.
%
%  F = normaltrace(B,d,w) computes the normal derivative of the field from 
%  discrete sources associated with bases B multiplied by weights w. Here B 
%  is an n vector of class basis and w is an n vector of scalar weights. 
%  F(i,j) is the field at point i from source j.
%
%  F = farfield(B,points,cof) computes the total far field at points given
%  by discrete sources associated with bases B, with the sources having
%  weights given by cof. Here B is an n vector of class basis and cof is an
%  m vector of scalar weights, where m is the total number of points
%  associated with all the bases in B. F(i) is the far field at point(i).
%
% See also: field, gradfield, farfield.
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


% Note: we make basis a child class of matlab.mixin.Heterogeneous so that
% we can apply the methods to arrays of basis.
classdef basis < handle & matlab.mixin.Heterogeneous
    
    properties
        kwave
        origin_sheet
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = basis(origin_sheet,kwave)
            
            % check that the origin_sheet is of the right class
            if ~isa(origin_sheet,'sheet_with_points')
                % ...if not display an error message
                error('origin_sheet must be of class sheet_with_points')
            end
            
            % store the inputs
            self.origin_sheet = origin_sheet;
            self.kwave = kwave;
            
        end
        
        %-----------------------------------------
        % compute the field from sources at the points
        % in origin_sheet at points on destination_sheet
        %-----------------------------------------
        
        function val = field(self,destination_sheet,weight)
            
            % Note: self is assumed to be an array
            
            % Specify default for the weights. Each entry in self should
            % have a weight.
            if nargin<3
                weight = ones(size(self));
            end
            
            %---------------------------------
            % setup
            %---------------------------------
            
            % Note: we assume that self is an array. First we need to get
            % the number of points on each associated sheet and put them in
            % n.
            % n(j) is the number of points in basis(j).origin_sheet
            
            % initialise n
            n = zeros(length(self),1);
            
            % store the number of points in each sheet
            for j=1:length(self)
                n(j) = size(self(j).origin_sheet.points,2);
            end
            
            % compute the cumulative sum... needed later
            sumn = [0;cumsum(n)];
            
            %---------------------------------
            % compute the field
            %---------------------------------
            
            % Initialise the field to the size of the destination sheet x
            % total number of source points. We will store
            % [field from basis(1),field from basis(2),...,field from
            % basis(end)]
            val = zeros(size(destination_sheet.points,2),sum(n));
            
            % loop through self
            for j=1:length(self)
                
                % in the appropriate part of val put the field from the
                % sources in basis(j)... multiplied by weight(j) and
                % multiplied by the quadrature weights
                val(:,sumn(j)+1:sumn(j+1)) = weight(j)...
                    * field(destination_sheet.points,...
                    self(j).origin_sheet.points,...
                    self(j).kwave) ...
                    * spdiags(self(j).origin_sheet.weights,0,size(self(j).origin_sheet.points,2),size(self(j).origin_sheet.points,2));
            end
            
        end
        
        %-----------------------------------------
        % compute the normal derivative of the field
        % from sources at the points in origin_sheet at
        % points on destination_sheet
        %-----------------------------------------
        
        function val = normaltrace(self,destination_sheet,weight)
            
            % Note: self is assumed to be an array
            
            % Specify default for the weights. Each entry in self should
            % have a weight.
            if nargin<3
                weight = ones(size(self));
            end
            
            %---------------------------------
            % setup
            %---------------------------------
            
            % Note: we assume that self is an array. First we need to get
            % the number of points on each associated sheet and put them in
            % n.
            % n(j) is the number of points in basis(j).origin_sheet
            
            % initialise n
            n = zeros(length(self),1);
            
            % store the number of points in each sheet
            for j=1:length(self)
                n(j) = size(self(j).origin_sheet.points,2);
            end
            
            % compute the cumulative sum... needed later
            sumn = [0;cumsum(n)];
            
            %---------------------------------
            % compute the field
            %---------------------------------
            
            % Initialise the field to the size of the destination sheet x
            % total number of source points. We will store
            % [field from basis(1),field from basis(2),...,field from
            % basis(end)]
            val = zeros(size(destination_sheet.points,2),sum(n));
            
            % compute the normal to the destination sheet at its mesh
            % points
            n = destination_sheet.sheet.normal(destination_sheet.pointsu,destination_sheet.pointsv);
            
            % loop through self
            for j=1:length(self)
                
                % compute the gradient of the field from the sources in 
                % basis(j) at the points on destination_sheet
                [dx,dy,dz] = gradfield(destination_sheet.points,...
                    self(j).origin_sheet.points,...
                    self(j).kwave);
                
                % put the quadrature weights into a diagonal sparse matrix
                A = spdiags(self(j).origin_sheet.weights,0,size(self(j).origin_sheet.points,2),size(self(j).origin_sheet.points,2));

                % apply the quadrature weights ny multiplying by the sparse
                % diagonal matrix
                dx = dx * A;
                dy = dy * A;
                dz = dz * A;
                
                % in the appropriate part of val put the field from the
                % sources in basis(j)... multiplied by weight(j).... do the
                % dot product with the normal by multiplying by each
                % component separately using sparse matrices.
                val(:,sumn(j)+1:sumn(j+1)) = weight(j)*(spdiags(n(1,:).',0,size(n,2),size(n,2)) * dx ...
                    + spdiags(n(2,:).',0,size(n,2),size(n,2)) * dy + spdiags(n(3,:).',0,size(n,2),size(n,2)) * dz);
                
            end
            
        end
        
        %-----------------------------------------
        % compute the total far field from weighted sources 
        % with weights for each source given by cof
        %-----------------------------------------
        
        function val = farfield(self,points,cof)
            
            % Note: self is assumed to be an array
            val = zeros(size(points,2),1);
            
            % loop through self
            for j=1:length(self)
                
                % add the far field from basis(j)
                val = val + farfield(points,...
                    self(j).origin_sheet.points,...
                    self(j).kwave) ...
                    * spdiags(self(j).origin_sheet.weights,0,size(self(j).origin_sheet.points,2),size(self(j).origin_sheet.points,2)) ...
                    * cof;
                
            end
            
        end
        
    end
    
end