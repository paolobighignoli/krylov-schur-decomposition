function [Us, Ts, complex] = sort_krylov_schur_decomposition(A, k)
%
%   move the k desired eigenvalues in the top left block of Ts
% 
%   INPUT
% 
%   A        matrix [n x n]
%   k        number of desired eigenvalues
%   cmd      comand to define the type of ordering
% 
%   OUTPUT
%
%   Us       ortonormal matrix [n x n]
%   Ts       matrix in Schur form [n x n]
%    
%   such that
%
%   Us' * A * Us = Ts 
%
%   where the desired eigenvalues are ordered on the top left block of Ts
%
%   complex = 0 if the last eigenvalue is real 
%   complex = 1 if the last eigenvalue is complex.
% 
    
    [U, T] = schur(A, 'real');
    
    autovT= ordeig(T);

    [~, i_autovordT] = sort(abs(autovT), 'descend');
   
    % questo interessa solo il caso abs e real 
    
    i1 = i_autovordT(k);% mi dice dove erano prima di essere ordinati
    i2 = i_autovordT(k+1);% cioè dove stanno in T 
    
    a1 = T(i1, i1);
    a2 = T(i2, i2);
    b1 = T(i2, i1);
    b2 = T(i1, i2);
    
    d = (a1 - a2)^2 + 4*b1*b2;
    
    if ((i2 - i1 == 1) && (d < 0)) % se k2 = k1+1 allora vuol dire che 
        complex = 1;% abbiamo una coppia complessa
    else 
        complex = 0;% nessuna coppia complessa
    end
    
   % nel caso complesso prende a coppie di default
   
    index = zeros(length(autovT), 1); % creo il vettore che mi dirà quali valori spostare
    index(i_autovordT(1:k+complex)) = true;% pone a 1 tutti i valori di interesse 
    
    [Us, Ts] = ordschur(U, T, index);
    
end