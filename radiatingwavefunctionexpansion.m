% Radiating wavefunction expansion
%
%  v = radiatingwavefunctionexpansion(u) copies a radiating wavefunction
%  expansion u into a new radiating wavefunction expansion v.
%
%  v = radiatingwavefunctionexpansion(u,x) creates a new radiating wavefunction
%  expansion v with origin x from the radiating wavefunction expansion u using the
%  translation addition theorem to implement the change of origin.
%
%  v = radiatingwavefunctionexpansion(u,x,n) creates a new radiating wavefunction
%  expansion v with order n and origin x from the radiating wavefunction 
%  expansion u using the translation addition theorem to implement the change 
%  of origin.
%
%  v = radiatingwavefunctionexpansion(n,x,k,cof) creates a wavefunction
%  expansion v with order n, expansion origin x, wavenumber k and
%  coefficient vector cof. The ordering of cof is
%
%    position    l     j
%    -------------------
%    1           0     0
%    2           1    -1
%    3           1     0
%    4           1     1
%    5           2    -2
%    6           2    -1
%    7           2     0
%    8           2     1
%    9           2     2
%
%  and so on, giving the coefficients for the wavefunction 
%  jh^{(1)}_l(k|x|) Y_{lj}(x/|x|}.
%
% Also:
%
%  f = u.evaluate(z) returns the values f of the wavefunction expansion u
%  evaluated at the points z, where z is a 3 x n array.
%
%  f = u.evaluate(z,mask) returns the values f of the wavefunction
%  expansion u evaluated at the points z for which mask==1 and NaN
%  elsewhere.
%
%  f = u.evaluateFarField(z) returns the values f of farfield of the
%  wavefunction expansion u evaluated at the points z on the unit circle, 
%  where z is a 3 x n array.
%
%  u.rotatecoordinatesabouty(theta) converts the wavefunction expansion to 
%  a new coordinate system obtained by rotation the old coordinate system
%  by angle theta about the y-axis. The expansion coefficients are
%  obtained accordingly.
%
%  u.rotatecoordinatesaboutz(theta) converts the wavefunction expansion to 
%  a new coordinate system obtained by rotation the old coordinate system
%  by angle theta about the z-axis. The expansion coefficients are
%  obtained accordingly.
%
%  c = u.getCoefficients() returns the vector of wavefunction expansion
%  coefficients.
%
%  u = v + w returns radiating wavefunction expansion obtained by adding
%  radiating wavefunction expansions v and w.
%
%  u = v - w returns radiating wavefunction expansion obtained by subtracting
%  radiating wavefunction expansion w from radiating wavefunction expansion v.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: wavefunctionexpansion, regularwavefunctionexpansion.
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


classdef radiatingwavefunctionexpansion < wavefunctionexpansion
    
    properties
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = radiatingwavefunctionexpansion(varargin)
            
            % check that we aren't being given an incident field... the
            % parent constructor would work with this but the coefficients
            % would not be correct.
            if nargin>=3 && isa(varargin{3},'incident')
                error('Cannot get radiating field coefficients from incident field')
            end
            
            % if we have one arguments then assume that varargin{1} is a
            % wavefunction expansion. Need to check that it is a radiating
            % wavefunction expansion
            if nargin==1 && ~isa(varargin{1},'radiatingwavefunctionexpansion')
                error('Given wavefunction expansion must be a radiatingwavefunctionexpansion')
            end
            
            % call parent constructor
            self = self@wavefunctionexpansion(varargin{:});

            % if varargin{1} is a wavefunction expansion and we have more 
            % arguments then assume that varargin{2} is a new
            % origin... so we need to modify the coefficients using the
            % translation addition theorem
            % If a third argument is provided then the order is going to
            % be changed too.
            if isa(varargin{1},'wavefunctionexpansion') && nargin > 1
                
                if isa(varargin{1},'radiatingwavefunctionexpansion')
                    
                    if nargin == 2
                        % only change the origin
                        self.apply_addition_theorem(varargin{2},'SAME');
                    else
                        % change the origin and the order
                        self.apply_addition_theorem(varargin{2},'SAME',varargin{3});
                        self.order = varargin{3};
                    end
                    
                else
                    
                    if nargin == 2
                        % only change the origin
                        self.apply_addition_theorem(varargin{2},'DIFF');
                    
                    else
                        % change the origin and the order
                        self.apply_addition_theorem(varargin{2},'DIFF',varargin{3});
                        self.order = varargin{3};                        
                    end
                end
                
            end
            
        end
    
        %-----------------------------------------
        % evaluate far field
        %-----------------------------------------
        
        function val = evaluateFarField(self,points)
            
            % get far field values
            val = sumcof(points,self.origin,self.kwave,self.coefficients,'F');
            
        end
        
        %-----------------------------------------
        % add
        %-----------------------------------------
        
        function val = plus(self,other)
            
            % check types match
            if ~isa(other,'radiatingwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = radiatingwavefunctionexpansion(self);
            
            % call parent add function
            plus@wavefunctionexpansion(val,other);
            
        end
        
        %-----------------------------------------
        % minus
        %-----------------------------------------
        
        function val = minus(self,other)
            
            % check types match
            if ~isa(other,'radiatingwavefunctionexpansion')
                error('Expansion types do not match')
            end
            
            % create a copy of self
            val = radiatingwavefunctionexpansion(self);
            
            % call parent add function
            minus@wavefunctionexpansion(val,other);
            
        end

    end
    
    methods(Access=protected)
        
        %-----------------------------------------
        % function to evaluate the expansion at given
        % points.
        % ** should be used only within methods
        % of the parent class. **
        %-----------------------------------------
        
        function val = internal_evaluate(self,points)
            
            val = sumcof(points,self.origin,self.kwave,self.coefficients,'H');
            
        end
        
    end
    
end