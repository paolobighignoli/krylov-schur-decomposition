function [E, W, is_complex, all_converged] = eigenvalues_krylov_schur(A, v1, n, k, m, maxIt, tol)
%
%  INPUT
%
%   A               matrix [n x n]
%   v1              initial vector [n X 1]
%   n               size of A
%   k               number of desired eigenvector
%   m               max dimension of Krylov decomposition
%   maxIt           max number of iteration
%   tol             tolerance
% 
%  OUTPUT
%
%   W               matrix of size [n x k] containing all the eigenvectors 
%   E               vector of size [n x 1] containing all the eigenvalues 
%   is_complex      = 0 if the last eigenvalue is real 
%                   = 1 if the last eigenvalue is complex
%   all_converged   all_converged = 0 if we reached expected number of converged eigenvalues
%                   all_converged = 1 if not
%

    [U, S, is_complex, all_converged, ~, ~, ~] = krylov_schur_decomposition(A, v1, k, m, maxIt, tol);
    [V, D] = eig(S(1:k+is_complex, 1:k+is_complex));
    
    E = diag(D);
    W = U(:, 1:k+is_complex) * V;
    
end

