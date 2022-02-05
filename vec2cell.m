% Convert vector to cell array
%
%  c = vec2cell(v) turns a vector of wavefunction expansion coefficients 
%  into a cell array.
%
%  This assumes that v holds coefficients v_{l,j} for l=0,...,n and j =
%  -l,...,l ordered thus:
%
%    index    l     j
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
%  and so on.
%
%  The cell array c stores the coefficients v_{l,-j},...,v_{l,j} in c{l+1}.
%
% See also: cell2vec.
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


function c = vec2cell(v)

% v should be of size (n+1)^2...
n = sqrt(size(v,1))-1;

% check that the size matches
if n~=round(n)
    error('v must be of length (n+1)^2 for some integer n')
end

% convert the cell
c = mat2cell(v,2*(0:n)+1,size(v,2));

