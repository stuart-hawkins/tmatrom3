% Wavefunction expansion
%
%  v = wavefunctionexpansion(u) copies a wavefunction expansion u into
%  a new wavefunction expansion v.
%
%  v = wavefunctionexpansion(n,x,inc) where inc is of type 'incident' 
%  creates a wavefunction expansion v with order n, expansion origin x and 
%  coefficients taken from inc.
%
%  v = wavefunctionexpansion(n,x,k,cof) creates a wavefunction expansion v 
%  with order n, expansion origin x, wavenumber k and coefficient vector 
%  cof.
%
% Also:
%
%  c = v.getCoefficients() returns the vector of wavefunction expansion
%  coefficients.
%
%  Note: this is an abstract class, it cannot be instantiated.
%
% See also: regularwavefunctionexpansion, radiatingwavefunctionexpansion.
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


classdef (Abstract) wavefunctionexpansion < tmatrom3
    
    properties
        order
        origin
        kwave
        coefficients
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = wavefunctionexpansion(varargin)
            
            if (nargin==1 || nargin==2 || nargin==3) && isa(varargin{1},'wavefunctionexpansion')
                                
                % v = wavefunctionexpansion(u)
                %
                % we are copying a given wavefunction expansion
                
                self.order = varargin{1}.order;
                self.origin = varargin{1}.origin;
                self.kwave = varargin{1}.kwave;
                self.coefficients = varargin{1}.coefficients(:);
                                
            elseif nargin==3
                
                % v = wavefunctionexpansion(n,x,inc)
                %
                % we are creating a wavefunction from an incident wave...
                % this *should* only be called by the regularwavefunction
                % expansion constructor
                
                self.order = varargin{1};
                self.origin = varargin{2};
                
                if ~isa(varargin{3},'incident')
                    
                    error('Third argument must be an incident wave')
                    
                end
                
                self.kwave = varargin{3}.kwave;
                
                self.coefficients = varargin{3}.get_coefficients(self.origin,...
                    self.order);
                
            elseif nargin==4
                
                % v = wavefunctionexpansion(n,x,k,cof)
                %
                % we are creating a wavefunction with the coefficients
                % given

                % check that the coefficient vector is of the correct size
                if length(varargin{4}) ~= (varargin{1}+1)^2
                    
                    error('The order n does not match the size of the coefficient vector cof.')
                    
                end
                
                self.order = varargin{1};
                self.origin = varargin{2};
                
                self.kwave = varargin{3};
                self.coefficients = varargin{4}(:);
                            
            else
                
                error('Incorrect parameters to constructor')
                
            end
            
        end
        
        %-----------------------------------------
        % evaluate near field
        %-----------------------------------------
        
        function val = evaluate(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % initialize array
            val = zeros(1,size(points,2));

            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
                        
            % apply mask if necessary
            if nargin>2
                points=points(:,mask);
            end
            
            % get far field values
            temp = self.internal_evaluate(points);
            
            % insert into return array
            if nargin>2
                val(mask) = temp;
                val(~mask) = NaN;
            else
                val = temp;
            end
                        
            % reshape the return array to match points
            if length(n)==2
                val = reshape(val,n(2),1);
            else
                val = reshape(val,n(2:end));
            end
            
        end
        
        %-----------------------------------------
        % get coefficients
        %-----------------------------------------

        function val = getCoefficients(self)
            
            val = self.coefficients;
            
        end

        %-----------------------------------------
        % rotate coordinates about z-axis
        %-----------------------------------------

        function varargout = rotatecoordinatesaboutz(self,angle)
            
            % convert the coefficients into a cell with c{l+1} holding the
            % coefficients of order l
            c = vec2cell(self.coefficients(:));

            for n=0:self.order

                % compute diagonal of rotation matrix
                alpha = exp(1i*(-n:n)*angle);

                % apply rotation matrix (which is diagonal)
                c{n+1} = alpha(:).* c{n+1};

            end

            % convert the coefficients from a cell to a vector
            self.coefficients = cell2vec(c);

            % return self if necessary
            if nargout > 0
                varargout{1} = self;
            end

        end
        
        %-----------------------------------------
        % rotate coordinates about y-axis
        %-----------------------------------------

        function varargout = rotatecoordinatesabouty(self,angle)

            % Wigner rotation matrices are represented by class wigner.
            % Setup an instanct of the wigner class
            wig = wavefunctionexpansion.get_wigner_matrix(self.order);

            % get the matrices that apply the rotation. W{l+1} applies the
            % rotation for harmonics of order l.
            W = wig.get(angle);

            % convert the coefficients into a cell with c{l+1} holding the
            % coefficients of order l
            c = vec2cell(self.coefficients(:));

            % loop through the orders
            for n=0:self.order

                % apply rotation matrix
                c{n+1} = W{n+1}.' * c{n+1};

            end

            % put the rotated coefficients back into a cell
            self.coefficients = cell2vec(c);

            % return self if necessary
            if nargout > 0
                varargout{1} = self;
            end

        end

        %-----------------------------------------
        % rotate coordinates about y-axis
        %-----------------------------------------

        function apply_addition_theorem(self,new_origin,type,new_order)
                        
            % set default for new_order
            if nargin<4
                new_order = self.order;
            end
            
            % get change of origin vector
            vec = new_origin - self.origin;
            
            % get polar coordinates for the translation vector
            r = sqrt(sum(vec.^2,1));
            theta = acos(vec(3)/r);
            phi = atan2(vec(2),vec(1));
            
            % We rotate the coordinate system so that the translation
            % vector is in the z-direction in the new coordinate system.
            % This reduces the complexity of applying the translation
            % addition theorem
            
            % rotate so that the translation vector is in the x-z plane
            % Note: this is changing the series coefficients according to
            % the rotation
            self.rotatecoordinatesaboutz(phi);
            
            % rotate again so that the translation vector is in the z
            % direction
            % Note: this is changing the series coefficients according to
            % the rotation
            self.rotatecoordinatesabouty(-theta);

            % get appropriate Bessel function
            if strcmp(type,'SAME')
                bess = sphbesselj(0:self.order+new_order,self.kwave*r);
            else
                bess = sphbesselh(0:self.order+new_order,self.kwave*r);
            end

            % in rotated coordinates theta = 0 so cos(theta) = 1
            Y = associatedLegendre(self.order+new_order,1);

            % get the wavefunction part of the translation addition theorem
            % coefficients
            for l=0:self.order+new_order
                c{l+1} = bess(l+1) * Y.get(l).';
            end

            % put the coefficients into a vector
            B = cell2vec(c(:));
            
            % apply the Wigner matrix to effect the transformation
            M = wavefunctionexpansion.get_transfer_matrix(new_order,self.order) * B;
            
            % update the coefficients
            self.coefficients = M * self.coefficients(:);
 
            % now we need to undo the rotation of the coordinate system
            % that we applied before applying the translation addition
            % theorem
            self.rotatecoordinatesabouty(theta);
            self.rotatecoordinatesaboutz(-phi);
  
            % update the origin
            self.origin = new_origin;
            
        end
        
        %-----------------------------------------
        % add
        %-----------------------------------------
        
        function plus(self,other)
            
            % check that the parameters match
            if self.order ~= other.order
                error('Orders do not match')
            end
            if max(abs(self.origin-other.origin)) > 1e-12
                error('Origins do not match')
            end
            if self.kwave ~= other.kwave
                error('Wavenumbers do not match')
            end
            
            % add the coefficients
            self.coefficients = self.coefficients + other.coefficients;
            
        end
        
        %-----------------------------------------
        % minus
        %-----------------------------------------
        
        function minus(self,other)
            
            % check that the parameters match
            if self.order ~= other.order
                error('Orders do not match')
            end
            if max(abs(self.origin-other.origin)) > 1e-12
                error('Origins do not match')
            end
            if self.kwave ~= other.kwave
                error('Wavenumbers do not match')
            end
            
            % add the coefficients
            self.coefficients = self.coefficients - other.coefficients;
            
        end

    end % end methods

    methods(Access=protected,Abstract=true)

        % this evaluates the wavefunction expansion at particular points...
        % it must be defined in the child class.
        val = internal_evaluate(self,points);
        
    end
    
    methods(Static)
        
        %-----------------------------------------
        % function to compute the transfer matrix
        %-----------------------------------------
        
        % This part is independent of the translation
        % and hence can be precomputed and stored. See Section 4 in
        % Dufva, Sarvas and Sten, Progress in Electromagnetics Research B,
        % Vol. 4, 79--99, 2008.
        %
        % Here we use a static method and persistent variables to store
        % the transfer matrix once and re-use it for all wavefunction
        % expansion objects. If a transfer matrix of higher order is
        % required then the higher order transfer matrix is stored. If a 
        % transfer matrix of lower order is subsequently required then the
        % appropriate subset is used. 
        % Note: this may lead to inconsistent accuracy but this is worth it 
        % for speed and reduced storage.
        
        function varargout = get_transfer_matrix(orderout,orderin)
        
            % declare persistent variables... these will behave a bit like 
            % static variables in other OO languages
            persistent nmax
            persistent transfer_matrix            

            % get order as larger of orderin and orderout
            order = max(orderin,orderout);
            
            % we use nmax being empty as a sign that the transfer matrix
            % has not been created yet.
            if isempty(nmax)
                
                if nargin < 1
                    error('order must be provided');
                else 
                    
                    % initialize nmax and transfer_matrix
                    nmax = 0;
                    transfer_matrix = [];
                
                end
                
            end
            
            % we assume we have the nmax transfer matrix. If the desired
            % order is bigger than nmax then we will need to recompute the
            % transfer matrix.
            if order ~= nmax
            
                % for brevity denote order by N
                N = order;
                
                % setup quadrature points
                [tp,pp,weights,points] = sphere_quadrature_points(2*N);

                % compute spherical harmonics at quadrature points                
                Y = associatedLegendre(2*N,cos(tp));
                
                % compute spherical harmonics at north pole                
                Y0 = associatedLegendre(2*N,1);

                % initialize matrix... matrix is actually a 3 dimensional
                % matrix and sparse. Matlab doesn't have 3D sparse matrices
                % so we define our own.
                transfer_matrix = sparse3((N+1)^2,(N+1)^2,(2*N+1)^2);
     
                % computation follows (46) in Dufva et al.
                for l=0:order
                    for j=-l:l

                        k = l^2 + l + j + 1;

                        for ld=abs(j):order

                            jd = j;

                            kd = ld^2 + ld + jd + 1;

                            for ldd = abs(l-ld):l+ld

                                jdd = jd-j;

                                kdd = ldd^2 + ldd + jdd + 1;

                                if abs(jdd)<=ldd && mod(l+ld+ldd,2)==0

                                    % use quadrature to evaluate the triple
                                    % integral
                                    int = sum(weights(:) ...
                                        .* Y.get(l,j) .* exp(-1i*j*pp(:)) ...
                                        .* Y.get(ld,jd) .* exp(1i*jd*pp(:)) ...
                                        .* Y.get(ldd,jdd) .* exp(-1i*jdd*pp(:)));

                                    transfer_matrix(k,kd,kdd) = (4*pi)*(-1i)^(ld-l-ldd)...
                                        * int;

                                end % if

                            end % ldd

                        end % ld

                    end % j

                end %l
            
                % store the size of the computed transfer matrix
                nmax = order;
                
            end

            % return the transfer matrix (or appropriate portion of it)
            if nargout > 0
                varargout{1} = transfer_matrix;
            end
            
        end
        
        %-----------------------------------------
        % function to setup the wigner matrix
        %-----------------------------------------

        % The Wigner matrix gives the coefficients for rotating the
        % coordinate system. We implement this with a class wigner that
        % wraps the matrix.
        
        function varargout = get_wigner_matrix(order)
        
            % declare persistent variables... these will behave a bit like 
            % static variables in other OO languages
            persistent wig
            
            if isempty(wig) || wig.order ~= order
            
                % initialise the wigner object
                wig = wigner(order);

                % setup the wigner object
                wig.setup();

            end

            if nargout > 0
                varargout{1} = wig;
            end

        end

    end

end