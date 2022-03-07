function [U, S] = reduce_krylov_schur_decomposition(U, S, k, m)
%
% reduce a Krylov decomposition of size m to one of size k (k < m)
% (you just keep the first k columns and last column of U and
%  the submatrix of dimension k of S plus his last row )
%
% INPUT
%
%   U             matrix of size [n x m+1]
%   S             matrix of size [m+1 x m]
%   k             final decomposition size
%   m             current decomposition size
%
% such that
%  
%           A * U(:, 1:m) = U * S 
%  
% OUTPUT
% 
%   U             matrix of size [n x k+1]
%   S             matrix of size [k+1 x k]
% 
% such that
%
%           A * U(:, 1:k) = U*S 
%
% or 
%
%           A * U(:, 1:k) = U(:, 1:k) * S(1:k, :) + U(:, k+1) * S(k+1, :)
%    

if nargin == 3
    m = size(S,2);
end

   U = [U(:, 1:k), U(:, m+1)];
   S = [S(1:k, 1:k); S(m+1, 1:k)];  
      
end                                     