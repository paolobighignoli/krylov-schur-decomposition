function flag = check_convergence(H, k, index, tol)
%
%   check if the eigenvalue in position index is converged,
%   as suggested in the paper we use the following criteria
%   
%   H(k+1,k) < max( norm(H) * epsilon, tol * lambda )
%   
%   INPUT
%   
%   H           matrix of size [k+1 x k] 
%   index       position of eigenvalue
%   k           size of decomposition
%   epsilon     machine precision 
%   tol         tolerance choosed by user 
%                   
%   OUTPUT
%
%   label        0 : not converged
%                1 : real eigenvalue converged
%                2 : complex eigenvalues converged
%
   
    epsilon = 2e-16;
    
    if index < k
        % calculate the discriminant to understand if I have a 
        % real or complex eigenvalue
        
        a1 = H(index,index);
        a2 = H(index+1,index+1);
        b1 = H(index+1,index);
        b2 = H(index, index+1);
        
        d = (a1-a2)^2 + 4*b1*b2; 
                                                        
    else 
        d = 1; 
    end
    
    if d > 0 % real case
        if abs(H(k+1, index)) < max(norm(H,'fro')*epsilon,abs(H(index, index))*tol )
            flag = 1;
        else 
            flag = 0;
        end
        
    else % complex case : delta < 0
        
        % calculate the complex eigenvalue
        lambda = (a1 + a2 + 1i*sqrt(-d)) / 2;
        
        if abs(H(k+1, index)) < max(norm(H, 'fro') * epsilon, abs(lambda) * tol )
            flag = 2;
        else
            flag = 0;
        end
    end
end