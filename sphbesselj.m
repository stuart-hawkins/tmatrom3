% Spherical Bessel function
%
%  J = sphbesselj(n,x) computes the spherical Bessel function of order n at
%  points x. Here n or x may be vectors.
%
%  See Abramowitz and Stegun, Handbook of Mathematical Functions, 
%  Section 10.1.1.
%
% See also: sphbesselh, sphbesseljd, sphbesselhd.
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


function val = sphbesselj(n,x)

val = sqrt(0.5*pi./x) .* besselj(n+0.5,x);
