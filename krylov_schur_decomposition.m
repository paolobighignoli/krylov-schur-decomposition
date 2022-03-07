function [U, S, is_complex, all_converged, n_converged_eigenvalues, n_iterations, H] = krylov_schur_decomposition(A, v_init, k, m, N, tol)
%
%   INPUT
%
%   A        matrix of size [n x n]
%   v_init   initial vector of size [n X 1]
%   n        size of A 
%   k        number of desired eigenvector
%   m        max dimension of Krylov decomposition
%   N        max number of iteration
%   tol      tolerance
% 
%  OUTPUT
%
%   U        orthogonal matrix of size [n x k+1]  
%            (or [n x k+2] if the last eigenvalue is complex)
%   S        Hessenberg matrix of dim. [k+1 x k] 
%            (or [k+2 x k+1] if the last eigenvalue is complex)
%
%   such that
% 
%   A * U(:, 1:k) = U * S 
%  (n, n) * (n, k) = (n, k+1) * (k+1, k)
% 
%   or equivalently
%
%   A * U(:, 1:k) = U(:, 1:k) * S(1:k, 1:k) + U(:, k+1) * S(k+1, :)
%  (n, n) * (n, k) = (n, k) * (k, k) + (n, 1) * (1, k)
%
%   is_complex                      = 0, if the last eigenvalue is real 
%                                   = 1, if the last eigenvalue is complex
%   all_converged                   = 0 if we reache expected number of converged eigenvalues
%                                   = 1 if not
%   n_converged_eigenvalues         number of converged eigenvalues
%   n_iterations                    number of iterations
%   H           history of all calculated eigenvalues
    
    n = size(A,1);
    
% set initial values
    U = zeros(n, m+1); 
    S = zeros(m+1, m); 
    H = zeros(k,1);
    is_complex = 0; 
    
% initialize first column of U
    U(:,1) = v_init / norm(v_init);
    
% create the starting Krylov decomposition from 0
    [U, S] = expand_krylov_schur_decomposition(A, U, S, 0, k); 
            
% do the loop until max number of iterations is reached or number of
% desired eigenvalues are calculated
    n_iterations = 0;
    p = 1;   
    while (n_iterations < N) && (p <= k) 
                              
        n_iterations = n_iterations+1;
        
        % expand krylov decomposition from k to m
        [U, S] = expand_krylov_schur_decomposition(A, U, S, k+is_complex, m);
        
        % order the desired eigenvalues 
        [W, T, is_complex] = sort_krylov_schur_decomposition(S(p:m, p:m), k-p+1);
        
        % save only the desired parts of our S, U matrices
        S(p:m, p:m) = T;
        S(1:p-1, p:m) = S(1:p-1, p:m) * W;
        U(:, p:m) = U(:, p:m) * W;
        S(m+1, p:m) = S(m+1, m) * W(end, :);
        
        % remove the part not desired
        [U, S] = reduce_krylov_schur_decomposition(U, S, k+is_complex, m); 
        
        % calculate the eigenvalues 
        H(1:k+is_complex,n_iterations) = eig(S(1:k+is_complex,1:k+is_complex));
        
        % test for convergence
        check = true; 
        while check
         
            % check if we have p converged eigenvalues
            result = check_convergence(S, k+is_complex, p, tol);
            
            if result  > 0
                %result == 1 || result == 2 
                % remember that result can be 0 (not converged), 1(real
                % converged), 2(complex pair converged)
                p = p + result; 
                
                if p > k 
                    check = false;
                end
            else
                check = false;
            end
        end
    end
    
    % return the convergence information
    if p > k
        all_converged = 0; % converged
        n_converged_eigenvalues = k + is_complex;
    else
        all_converged = 1; % not converged
        n_converged_eigenvalues = p - 1;
    end
    
    H(1:k+is_complex,n_iterations) = eig(S(1:k+is_complex,1:k+is_complex));
        
end