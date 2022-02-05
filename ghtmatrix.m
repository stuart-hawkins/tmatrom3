% Create the T-matrix from far-field data.
% 
%  T = tmatrix(n,k,s) computes the T-matrix of order n for wavenumber k
%  and origin [0;0;0]. The scatter is completely described by the object s
%  which must be of class solver.
%
%  T = tmatrix(n,k,s,x0) returns the T-matrix with origin specified by x0.
%
% This code uses the numerically stable algorithm in Ganesh and Hawkins,
% ANZIAM J. (50), p. C121--C136, 2008.
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


function tmat = ghtmatrix(order,kwave,solver,origin)

% set default for origin
if nargin < 4
    origin = [0;0;0];
end

% check that solver is of the correct type
if ~isa(solver,'solver')
    error('solver must be an instance of the solver class')
end

% - - - - - - - - - - - - - - - - -
% setup the quadrature points
% - - - - - - - - - - - - - - - - -

[tp,pp,weight,points] = sphere_quadrature_points(order);

% - - - - - - - - - - - - - - - - -
% solve the scattering problems
% - - - - - - - - - - - - - - - - -

% setup a cell array of incident fields
kk = 1;
for n = 0:order
    for j = -n:n
        inc{kk} = regularwavefunction(n,j,kwave,origin);
        kk = kk + 1;
    end
end

% set the incident field in the solver
solver.setIncidentField(inc);

% solve the scattering problems
solver.solve();

% get the farfield
farfield = solver.getFarField(points);

% assume the farfield was computed with origin 0... we need
% to adjust for the modified origin
if max(abs(origin)) > 0
    
    error('not yet implemented')
    
end


% - - - - - - - - - - - - - - - - -
% compute the T-matrix
% - - - - - -- - - - - - - - - - -

% Note: we should really split this into two parts similar to in Ganesh and
% Graham, Journal of Computational Physics, 198 p. 211-242 2004.
% However the CPU time for this part is small so this direct approach
% suffices.

% get the associate Legendre values
Y = associatedLegendre(order,cos(tp(:)));

% get the spherical harmonic in a cell array... tmp{n+1} contains the
% values of the spherical harmonic of order n
for n = 0:order
    
    tmp{n+1} = (Y.get(n) .* exp(-1i*pp(:)*(-n:n))).';
    
end

% convert the cell array into a vector
H = cell2vec(tmp(:));

% get a matrix containing the values of n... this is used for the scaling
% in the next part
for n=0:order
    nc{n+1} = n*ones(2*n+1,1);
end
nn = cell2vec(nc(:));

% compute the T-matrix using (16) in Ganesh and Hawkins
% ANZIAM J. 51 C215--C230 (2010)
matrix = diag(kwave*(1i).^(nn+1)) * H * diag(weight) * farfield;
    
% put some basic information in the comment string
comment = sprintf('Computed using ghtmatrix v%0.1f with solver %s. %s',tmatrom3.version(),class(solver),solver.description());

% creat the T-matrix object
tmat = tmatrix(order,kwave,matrix,origin,comment);
