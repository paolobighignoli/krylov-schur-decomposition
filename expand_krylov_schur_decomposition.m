function [U, S] = expand_krylov_schur_decomposition(A, U, S, k, m)
% 
% expand a Krylov decomposition of level k to one of level m
% where k<m
%
% INPUT
%
% A matrix [n x n]
% Q matrix [n x k+1]
% H matrix [k+1 x k]
% k current decomposition size
% m final decomposition size
%
% such that
%
% A * Q(:,1:k) = Q(:,1:k+1) * H(1:k+1,1:k)
% A * Q(:,1:k) = Q(:,1:k) * H(1:k,1:k) + Q(:,k+1)*H(k+1,:)
%
% OUTPUT
% 
% Q matrix [n x m+1]
% H matrix [m+1 x m]
% 
% such that
%
% A * Q(:,1:m) = Q(:,1:m+1) * H(1:m+1,1:m)
% A * Q(:,1:m) = Q(:,1:m) * H(1:m,1:m) + Q(:,m+1)*H(m+1,:)
%

    for j = k+1 : m
        v = A*(U(:, j));
        
        s1 = U(:, 1:j)' * v; % first normalization
        v = v - U(:, 1:j) * s1;
        
        s2 = U(:, 1:j)' * v; % second normalization
        v = v - U(:, 1:j) * s2;
        
        s = s1 + s2;
        
        U(:, j+1) = v / norm(v);       
        S(1:j+1, j) = [s; norm(v)];
    end

end

