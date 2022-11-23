function [C] = corr_col(A, B, dim)
%CORR_COL quickly compute the column by column correlation between two large matrices
%
% C=CORR_COL(A, B) computes the correlation for each column in A of the
% congruent column in matrix B.  A and B must be identical in size
%
% C=CORR_COL(A, B, dim) computes the correaltions for the dimension
% specified in dim
%
% This function produces nearly identical values as diag( corr( A,B ))
% but can be run on matrices that might be too large for diag( corr( A,B ))
%
% See also diag corr
%
%   The most up to date version of this code can be found at:
%   https://github.com/slayton/matlab-corr-col
%
% Copyright(c) 2012, Stuart P. Layton <stuart.layton@gmail.com> MIT
% http://stuartlayton.com



% check to see if A and B are both matrices
if ~ismatrix(A) || ~ismatrix(B)
   error('Input matrices A,B can only have 2 dimensions'); 
end

% Ensure that A and B are the same size
if ~all( size(A) == size(B) )
    error('A and B must be the same size');
end

%If dim isn't specified default 1
if nargin == 2
    dim = 1;
end

if dim > ndims(A) || dim<1
    error('Invalid dimension specified, dim must be between 1 and the number of dimensions on A');
end


% Make each column zero-mean
A = bsxfun( @minus, A, mean( A, dim) );
B = bsxfun( @minus, B, mean( B, dim) );

% L2 normalize each column
A = bsxfun( @times, A, 1./sqrt( sum( A.^2, dim) ) );
B = bsxfun( @times, B, 1./sqrt( sum( B.^2, dim) ) );

% Take the dot product of the columns and then sum
C=sum( A.*B, dim);

