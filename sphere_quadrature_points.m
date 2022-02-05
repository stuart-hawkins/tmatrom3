% Gauss-Legendre quadrature points on the sphere
%
%  [theta,phi] = sphere_quadrature_points(n) returns the spherical polar
%  coordinates (theta,phi) for 2*(n+1)^2 Gauss-Legendre quadrature points on 
%  the surface of the unit sphere.
%
%  [theta,phi,weights,points] = sphere_quadrature_points(n) additionally
%  returns a 3 x 2*(n+1)^2 matrix points containing the Cartersian 
%  coordinates for the Gauss-Legendre points, and a vector weights
%  containing the corresponding quadrature weights.
%
%  See Section 3 of Ganesh and Hawkins, Numerical Algorithms 43:25-60
%  (2006) for details of the Gauss-Legendre rule.
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


function [theta,phi,weights,points] = sphere_quadrature_points(order)

% calculate the phi points (azimuth angle)
pp = pi*(0:2*order+1).'/(order+1);
pw = pi*ones(2*order+2,1)/(order+1);

% calculate the theta points (polar angle)
[tp,tw]=gauss_legendre(order);

% gauss_legendre does probabilists version... multiply by 2 for
% physicists version
tw = 2*tw;

% ensure these are column vector
tp = acos(tp(:));

% get the lengths of the arrays
n = length(tp);
m = length(pp);

% make a grid of the theta points and weights
tp = repmat(tp(:).',m,1);
tw = repmat(tw(:).',m,1);

% make a grid of the phi points and weights
pp = repmat(pp(:),1,n);
pw = repmat(pw(:),1,n);

% the final weights are the product of the theta and phi weights
ww = tw.*pw;

% turn the matrices of points and weights into column vectors
theta = tp(:).';
phi = pp(:).';
weights = ww(:).';

% compute the corresponding Cartesian coordinates
points = [sin(theta).*cos(phi);sin(theta).*sin(phi);cos(theta)];
