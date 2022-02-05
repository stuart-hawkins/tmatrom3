% Farfield of the Greens function for the Helmholtz equation
%
%  v = farfield(p,q,k) computes the farfield of the Green's function for 
%  the Helmholtz equation with wavenumber k for sources q at points p. If 
%  p is a 3 x n array and q is a 3 x m array then v is an n x m array with 
%  v(i,j) being the field at p(:,i) from a source at q(:,j).
%
% See also: field, gradfield.
%
% Stuart C. Hawkins - 22 September 2021

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


function val = farfield(p,q,kwave)

% get the sizes of p and q
n = size(p,2);
m = size(q,2);

% px is an n x m array with px(i,j) = p(1,i)
px = repmat(p(1,:).',1,m);

% py is an n x m array with px(i,j) = p(2,i)
py = repmat(p(2,:).',1,m);

% pz is an n x m array with px(i,j) = p(3,i)
pz = repmat(p(3,:).',1,m);

% qx is an n x m array with qx(i,j) = q(1,j)
qx = repmat(q(1,:),n,1);

% qy is an n x m array with qx(i,j) = q(2,j)
qy = repmat(q(2,:),n,1);

% qz is an n x m array with qx(i,j) = q(3,j)
qz = repmat(q(3,:),n,1);

% evaluate the dot product p.q
dp = px.*qx + py.*qy + pz.*qz;

% compute the farfield of the Green's function
val = (0.25/pi)*exp(-1i*kwave*dp);
