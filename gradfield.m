% Gradient of the Greens function for the Helmholtz equation
%
%  [dx,dy,dz] = gradfield(p,q,k) computes the gradient of the Green's function 
%  for the Helmholtz equation with wavenumber k for sources q at points p. 
%  If p is a 3 x n array and q is a 3 x m array then dx is an n x m array 
%  with dx(i,j) being the x-component of the gradient at p(:,i) from a 
%  source at q(:,j). Similarly dy and dz contain the y and z coordinates of
%  the gradient respectively.
%
% See also: field.
%
% Stuart C. Hawkins - 23 August 2021

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


function [dx,dy,dz] = gradfield(p,q,kwave)

% get the sizes of p and q
n = size(p,2);
m = size(q,2);

% dx is an n x m array with dx(i,j) = p(1,i)-q(1,j)
dx = repmat(p(1,:).',1,m)-repmat(q(1,:),n,1);

% dy is an n x m array with dy(i,j) = p(2,i)-q(2,j)
dy = repmat(p(2,:).',1,m)-repmat(q(2,:),n,1);

% dz is an n x m array with dz(i,j) = p(3,i)-q(3,j)
dz = repmat(p(3,:).',1,m)-repmat(q(3,:),n,1);

% r is an n x m array with r(ij) = |p(:,i)-q(:,j)|.
r = sqrt(dx.^2 + dy.^2 + dz.^2);

% evaluate the scalar part... d/dr of the Green's function x 1/r
tmp = (1i*kwave*r-1).*exp(1i*kwave*r)./(4*pi*r.^3);

% get the vector components by multiplying tmp with vector r
dx = tmp.*dx;
dy = tmp.*dy;
dz = tmp.*dz;

