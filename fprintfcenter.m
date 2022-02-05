% Display centred text
%
%   fprintfcenter(fmt,a,...) displays a with the format fmt similar to the 
%   fprintf function, but centred in the command window.
% 
% See also: fprintf.
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


function val = fprintfcenter(fmt,varargin)

% work out the width in characters of the command window
tmp = get(0,'CommandWindowSize');
width = tmp(1);

% format the output using the sprintf command
text = sprintf(fmt,varargin{:});

% get the length of the formatted text
len = length(text);

% work out the padding that is needed (in characters) on the left to center
% the text
offset = floor((width-len)/2);

% display the text... apply the padding as a string generated using repmat
val = fprintf('%s%s\n',repmat(' ',offset,1),text);
