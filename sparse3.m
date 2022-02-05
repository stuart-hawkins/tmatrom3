% Sparse 3D matrix
%
%  M = sparse(n,m,p) initialises an n x m x p sparse matrix.
%
%  M(i,j,k) returns the (i,j,k) entry of M.
% 
%  M(i,j,k) = x sets the (i,j,k) entry of M to x.
%
%  A = M * x computes the matrix-vector product of M with x. Here the
%  multiplication is with respect to the third dimension of M and the
%  result is an n x m matrix.
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


classdef sparse3 < handle
    
    properties
        i
        j
        k
        entry
        n
        m
        p
    end

    methods

        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = sparse3(n,m,p)

            % store size
           self.n = n;
           self.m = m;
           self.p = p;
           
           % ...other properties are empty, signifying the matrix has no
           % entries at the moment
           self.i = [];
           self.j = [];
           self.k = [];
           self.entry = [];
           
           % Note: the (self.i(s),self.j(s),self.k(s)) entry of the matrix
           % is self.entry(s).

        end

        %-----------------------------------------
        % overload subscript assignment
        %-----------------------------------------

        function out = subsasgn(self,S,val)

            % create a new matrix for the output
            out = sparse3(self.n,self.m,self.p);

            % check that we are doing a self(i,j,k) type assignment, and
            % that (i,j,k) are within the dimensions of the matrix
            if strcmp(S.type,'()') && S.subs{1} <= self.n ...
                    && S.subs{2} <= self.m ...
                    && S.subs{3} <= self.p 
                
                % ... store the entry

                % Note: we don't mind duplicated (i,j,k)s because we only
                % want this matrix for matrix-vector multiplication, which
                % will come out right for duplicated (i,j,k)s
                
                % add the new entry to the (i,j,k) arrays
                out.i = [self.i;S.subs{1}];
                out.j = [self.j;S.subs{2}];
                out.k = [self.k;S.subs{3}];
                
                % add the new value 
                out.entry = [self.entry;val];

            else

                % ... don't want to do anything - this kind of indexing is
                % not supported. Just copy the matrix.
                out.i = self.i;
                out.j = self.j;
                out.k = self.k;
                out.entry = self.entry;

            end

        end

        %-----------------------------------------
        % overload subscript referencing
        %-----------------------------------------

        function val = subsref(self,S)

            % check that we are doing S(i,j,k) type referencing
            if strcmp(S.type,'()')

                % find the entries matching the given (i,j,k)
                if length(self.i) > 0
                    kk = find(self.i==S.subs{1} & self.j==S.subs{2} & self.k==S.subs{3});
                else
                    kk = [];
                end
                
                % Note: there may be many matching entries because we allow
                % duplicate (i,j,k)s.... just add up duplicates.
                if ~isempty(kk)
                    val = sum(self.entry(kk));
                else
                    val = 0;
                end

            end

        end

        %-----------------------------------------
        % matrix vector product
        %-----------------------------------------

        function val = mtimes(self,x)

            % check that the size of the matrix and the size of x match
            if self.p~=size(x,1)
                error('sparse3 inner dimension mismatch')
            end

            % check that x is a column vector
            if size(x,2)~=1
                error('x must be a vector')
            end

            % initialise the result.... note the result will be an n x m
            % matrix.
            val = zeros(self.n,self.m);

            % loop through entries in the matrix
            for kk=1:length(self.i)
                % do the multiplication
                val(self.i(kk),self.j(kk)) = val(self.i(kk),self.j(kk)) ...
                    + self.entry(kk) * x(self.k(kk));
            end

        end

    end

end