% TMATROM base class
%
%  This class is intended as a base class for the main TMATROM classes,
%  such as tmatrix and wavefunctionexpansion etc. It handles version
%  numbering and the splash message.
%
%  tmatrom3.version() returns the version number.
%
%  tmatrom3.date() returns the date of the current version
%
% See also: tmatrix, wavefunctionexpansion.
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


classdef tmatrom3 < handle
    
    properties        
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = tmatrom3()
           
            % call initialise message... this displays the version number
            % etc
            tmatrom3.initialise();
            
        end
        
    end
    
    methods(Static,Sealed)
        
        %-----------------------------------------
        % version number for package
        %-----------------------------------------

        function val = version()
            
            val = 1.0;
            
        end
        
        %-----------------------------------------
        % date for current release of the package
        %-----------------------------------------

        function val = date()
            
            val = '22 September 2021';
            
        end

        %-----------------------------------------
        % print version and copyright
        % information for the package
        %-----------------------------------------
        
        function showInfo()
            
            fprintfcenter('TMATROM3 (acoustic)');
            fprintfcenter('Version %0.1f (%s)',tmatrom3.version(),...
                tmatrom3.date());
            fprintfcenter('%s','(c) Stuart C. Hawkins and M. Ganesh, 2021');
            
        end
        
        %-----------------------------------------
        % print status information for the package.
        %-----------------------------------------
        
        function showStatus()

            % nothing to do... no things saved in current version
            
        end
        
        %-----------------------------------------
        % acts like a static variable that stores
        % whether or not the class has been used.
        % Mainly this is used to show the info message
        % the first time the class is used.
        %-----------------------------------------
        
        function initialise(varargin)
            
            persistent has_been_initialised
            
            % if class has not been initialised yet...
            if isempty(has_been_initialised) || has_been_initialised==0
                
                if nargin==0
                    
                    % Note: if a parameter is provided then the status
                    % message is not printed. Thus the welcome message can
                    % suppressed if desired.
                    
                    % ...then show class info and status
                    tmatrom3.showInfo()
                    fprintf('\n');
                    tmatrom3.showStatus()
                    
                end
                
                % record that the class has been initialised
                has_been_initialised = 1;
                
            end
            
        end
        
    end % static methods
    
end
        