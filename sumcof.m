% Evaluate wavefunction expansion series.
%
%   z = sumcof(x,x0,k,c,'H') returns the values z of the radiating 
%   wavefunction expansion with coefficients c, centre x0 and wavenumber k 
%   at points x. Here x should be a 3 x n matrix.
%
%   z = sumcof(x,x0,k,c,'J') returns the values z of the regular 
%   wavefunction expansion with coefficients c, centre x0 and wavenumber k 
%   at points x. Here x should be a 3 x n matrix.
%
%   z = sumcof(x,x0,k,c,'F') returns the values z of the far field of the 
%   radiating wavefunction expansion with coefficients c, centre x0 and 
%   wavenumber k at points x on the unit sphere. Here x should be a 
%   3 x n matrix.
%
% Note: in the above vectors in the plane are represented by
% complex numbers.
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


function val = sumcof(points,centre,kwave,cof,type)

%-------------------------------------------------
% setup
%-------------------------------------------------

% make sure the coefficient vector is a column vector
cof=cof(:);

% make sure the centre is a column vector
centre = centre(:);

% determine the maximum order from the length of the coefficient vector
nmax = sqrt(length(cof))-1;

% report an error if the coefficients vector isn't the right size
if nmax~=round(nmax)
    error('cof must be of length (n+1)^2 for some integer n')
end

% create a vector of indexes.... helps to vectorize the computation
n=0:nmax;

%-------------------------------------------------
% turn points into a vector
%-------------------------------------------------

% get the size of points... we need this later to reshape the result
mp = size(points);

% reshape points into a 3 x [] array
points = reshape(points,3,prod(mp(2:end)));

% get the shape of points so we can restore it later
np=size(points);

% check dimensions of points are okay
if np(1)~=3 || length(np)~=2
    error('points must be a 3 x n array')
end

% subtract the local origin if required
if strcmp(type,'F')
    % Note: for far field the far field is computed with respect to the
    % true origin
    p = points;
else   
    p = points - repmat(centre,1,np(2));
end

%-------------------------------------------------
% compute the field
%-------------------------------------------------

% convert to polar coordinates
rad = sqrt( p(1,:).^2 + p(2,:).^2 + p(3,:).^2 );
theta = acos(p(3,:)./rad);
phi = atan2(p(2,:),p(1,:));

% make rad, theta and phi column vectors
rad = rad(:);
theta = theta(:);
phi = phi(:);

% make a matrix from n and rad
[nd,rd]=meshgrid(n,kwave*rad);

% get Bessel/Hankel/far field values as appropriate
if strcmp(type,'J')
    
    bess=sphbesselj(nd,rd);
    
elseif strcmp(type,'H')
    
    bess=sphbesselh(nd,rd);
    
elseif strcmp(type,'F')
    
    bess=(1/kwave)*(-1i).^(nd+1);   
    
end

% compute the angular part
Y = associatedLegendre(nmax,cos(theta));

% adjust the farfield if the centre of the scatterer is not
% the origin... there needs to be a change of phase
if strcmp(type,'F') && max(abs(centre))~=0
    dp = centre.' * p;
    phase = exp(-1i*kwave*dp);
    bess = diag(phase)*bess;
end

% get the coefficients in the form of a cell array.... ccof{l+1} holds the
% coefficients of order l
ccof = vec2cell(cof);

% initialise output
val = zeros(np(2),1);

% loop through the orders
for n=0:nmax
    % compute the sum of order n terms
    val = val + bess(:,n+1) .* ((Y.get(n) .* exp(1i*phi*(-n:n))) * ccof{n+1});
end

%-------------------------------------------------
% make the return value the same shape as the original
% array of points
%-------------------------------------------------

% compute the sum of the wavefunctions and reshape
if length(mp)==2
    val=reshape(val,mp(2),1);
else
    val=reshape(val,mp(2:end));
end