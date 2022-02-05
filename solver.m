% Solver
%
%  s = solver(k,inc) creates a solver object with wavenumber k and incident
%  fields specified in the incident field inc.
%
%  s = solver(k,[]) creates a solver object with wavenumber k and defers
%  setting the incident fields.
%
%  s.setIncidentField(inc) sets the incident field as specified in the cell
%  array inc.
%
%  val = s.getCrossSection(pts) computes the cross section for the
%  incident field inc{1} at points on the unit sphere specified in a 3 x n 
%  matrix pts.
%
%  val = s.getCrossSection(pts,index) computes the cross section for the
%  incident fields inc{index} at points on the unit sphere specified in a 
%  3 x n matrix pts.
%
% Note: this is an abstract class that provides common methods for solvers
% and specifies the interfaces for other methods that must be provided in
% the child class.
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


classdef solver < handle
    
    properties
        kwave
        incidentField
        numIncidentField
    end
    
    methods
        
        %===============================================================
        % methods that will (probably) not be overridden in child class
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solver(kwave,incidentField)
            
            self.kwave = kwave;
            
            if ~iscell(incidentField)
                
                % turn a 'scalar' incident field into a length 1 cell
                self.incidentField{1} = incidentField;
               
            else
                
                % incident field is already given as a cell
                self.incidentField = incidentField;
                
            end
            
            % set the number of incident fields
            self.numIncidentField = length(self.incidentField);
            
        end
        
        %-----------------------------------------
        % set incident field
        %-----------------------------------------
        
        function setIncidentField(self,incidentField)
            
            % make sure that the incident field is stored as a cell array
            % even if there is only one incident field
            if ~iscell(incidentField)
                
                % clear the incident field
                self.incidentField = [];
                
                % turn a 'scalar' incident field into a length 1 cell                
                self.incidentField{1} = incidentField;
               
            else
                
                % incident field is already given as a cell
                self.incidentField = incidentField;
                
            end
            
            % set number of incident fields
            self.numIncidentField = length(self.incidentField);
            
        end
        
        %-----------------------------------------
        % get cross section in dB
        %-----------------------------------------

        function val = getCrossSection(self,points,index)
            
            % set default for index
            if nargin < 3
                index = 1;
            end
            
            val = 10*log10(4*pi*abs(self.getFarField(points,index)).^2);
            
        end
        
    end
            
    methods(Abstract=true)
        
        %===============================================================
        % these methods must be overridden in the child class
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        setup(self)
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        solve(self)
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        val = getFarField(self,points,index)
        
        %-----------------------------------------
        % return a string describing the solver
        %-----------------------------------------
        
        val = description(self)

    end % end methods
    
end