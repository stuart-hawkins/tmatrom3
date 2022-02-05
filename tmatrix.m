% T-matrix class
%
%   T = tmatrix(fname) loads a T-matrix from file fname. The file is
%   assumed to be of Matlab *.mat format.
%
%   T = tmatrix(n,k,M) creates a T-matrix object of order n, wavenumber k, 
%   and matrix M. The origin is set to be the default [0;0;0].
%
%   T = tmatrix(n,k,M,x) creates a T-matrix object of order n, wavenumber k, 
%   matrix M and origin x.
%
%   T = tmatrix(n,k,M,x,str) creates a T-matrix object of order n, wavenumber k, 
%   matrix M, origin x and comment str. 
%
%   S = tmatrix(T) copies T into a new tmatrix object.
%
% Also:
%
%   val = T.error() gives a measure of the error in the T-matrix based on a
%   symmetry relation. See Equation (7.67) in P. A. Martin, Multiple 
%   Scattering: Interaction of Time-Harmonic Waves with N obstacles for
%   details. This symmetry relation is valid for non-absorbing obstacles
%   and may not be correct for absorbing particles (complex refractive
%   index, or Robin boundary conditions).
% 
%   T.setOrigin(x) sets the origin of the T-matrix to x. This is a virtual
%   origin that is only used in interactions with wave functions. This does
%   not change the T-matrix.
%
%   T.setComments(str) stores str as comments associated with the T-matrix.
%
%   T.getComments() prints any comments associated with the T-matrix.
%
%   val = T.getComments() returns any comments associated with the T-matrix
%   in val.
%
%   T.save(fname) saves the T-matrix in .mat format.
%
% Example:
%
%   p = plane_wave(0,k);
%   u = regularwavefunctionexpansion(n,0,p);
%   T = tmatrix('sample.mat');
%   v = T * u;
%
%   Now v is a radiating wavefunction expansion for the scattered field
%   induced by the plane wave p. The T-matrix is loaded from file
%   sample.mat.
%
% See also: regularwavefunctionexpansion, radiatingwavefunctionexpansion,
% plane_wave, point_source, ghtmatrix.
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


classdef tmatrix < tmatrom3
    
    properties
        order
        kwave
        origin
        matrix
        comments
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = tmatrix(varargin)
            
            if nargin==1
                
                if isa(varargin{1},'tmatrix')
                    
                    % then we are copying a T-matrix
                    self.order = varargin{1}.order;
                    self.kwave = varargin{1}.kwave;
                    self.matrix = varargin{1}.matrix;
                    self.origin = varargin{1}.origin;
                    self.comments = varargin{1}.comments;
                
                else
                    
                    % we are loading a T-matrix
                    self.load(varargin{1})
                    
                end
                    
            elseif nargin==2
                
                % then we are copying a T-matrix
                % but changing the origin
                self.order = varargin{1}.order;
                self.kwave = varargin{1}.kwave;
                self.matrix = varargin{1}.matrix;
                self.origin = varargin{2};
                self.comments = varargin{1}.comments;
                
            elseif nargin==3
            
                % then we are getting the T-matrix from a matrix with the
                % default origin
                self.order = varargin{1};
                self.kwave = varargin{2};
                self.matrix = varargin{3};
                self.origin = [0;0;0];
            
            elseif nargin==4
                
                % then we are getting the T-matrix from a matrix
                self.order = varargin{1};
                self.kwave = varargin{2};
                self.matrix = varargin{3};
                self.origin = varargin{4};

            elseif nargin==5
                
                % then we are getting the T-matrix from a matrix and a
                % comment is supplied
                self.order = varargin{1};
                self.kwave = varargin{2};
                self.matrix = varargin{3};
                self.origin = varargin{4};
                self.comments = varargin{5};
                
            end

        end
        
        %-----------------------------------------
        % multiply tmatrix x regular wave expansion
        %-----------------------------------------

        function val = mtimes(self,expansion)
            
            % - - - - - - - - - - - - - - - - - 
            % check the T-matrix and the wave 
            % expansion are compatible
            % - - - - - - - - - - - - - - - - - 

            if ~isa(expansion,'regularwavefunctionexpansion')
                
                error('expansion must be a regularwavefunctionexpansion')
                
            end
            
            if self.kwave ~= expansion.kwave
                
                error('T-matrix and expansion wavenumbers do not match.')
                
            end
            
            if max(abs(self.origin-expansion.origin)) > 0
                
                error('T-matrix and expansion centers do not match.')
                
            end
           
            if self.order ~= expansion.order
                
                error('T-matrix and expansion orders do not match.')
                
            end
            
            % - - - - - - - - - - - - - - - - - 
            % do product 
            % - - - - - - - - - - - - - - - - - 

            % create a radiating wave function expansion with coefficients
            % obtained by matrix multiplication with the T-matrix
            val = radiatingwavefunctionexpansion(self.order,self.origin,...
                self.kwave,self.matrix * expansion.coefficients(:));
            
        end
        
        %-----------------------------------------
        % error check
        %-----------------------------------------

        % Check the symmetry relation Equation (7.67) in P. A. Martin, 
        % Multiple Scattering: Interaction of Time-Harmonic Waves with N 
        % obstacles for details.

        function val = error(self,opts)
           
            val = max(max(abs(self.matrix + self.matrix' ...
                + 2 * self.matrix' * self.matrix)));
            
            if nargin>1
                val = val / max(max(abs(self.matrix)));
            end
            
        end

        %-----------------------------------------
        % load
        %-----------------------------------------

        function load(self,fname)
            
            % put the class properties into a struct
            data = load(fname);
            
            % ** in future released we might need to check the version here
            % eg if class properties change **
            
            % set the class properties from the struct
            self.order = data.order;
            self.kwave = data.kwave;
            self.origin = data.origin;
            self.matrix = data.matrix;
            self.comments = data.comments;
            
        end
        
        %-----------------------------------------
        % save
        %-----------------------------------------

        function save(self,fname)
            
            % put the class properties in a struct
            data = struct('order',self.order,'kwave',self.kwave,...
                'origin',self.origin,'matrix',self.matrix,...
                'version',self.version(),'comments',self.comments);
            
            % save the struct
            save(fname,'-struct','data');
            
        end
        
    end
    
end