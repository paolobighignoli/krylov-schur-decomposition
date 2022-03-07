

% final test 

rng(123);

N = 100;
m = 30;
k = 20;
tol = 2e-16;

% A = rand(N,N);
% A = gallery('frank', N)*100;
% A = gallery('circul', N);
% A = gallery('cycol', N);
% A = gallery('cycol', N);
A = magic(100);
v = rand(N,1);



all_eigenvalues = eig(A);
[~, index] = sort(abs(all_eigenvalues),'descend');
sorted_eigenvalues = all_eigenvalues(index);
sorted_k_eigenvalues = sorted_eigenvalues(1:k);
% sorted_k_eigenvalues = sorted_eingevalues(1:k);


% plot(real(all_eigenvalues),imag(all_eigenvalues),'b*',real(sorted_k_eigenvalues),imag(sorted_k_eigenvalues),'ro');

%% PLOT KRYLOV-SCHUR

[eigenvalues_krylov, eigenvectors_krylov, isC, flag] = eigenvalues_krylov_schur(A, v, N, k, m, 2000, 2e-16);

plot(real(all_eigenvalues),imag(all_eigenvalues),'b*',...
real(eigenvalues_krylov),imag(eigenvalues_krylov),'ro',...
real(sorted_k_eigenvalues),imag(sorted_k_eigenvalues),'g*');

%% PLOT EIGENVALUES IN MOTION

[~, ~, ~, ~, ~, ~, H] = krylov_schur_decomposition(A, v, k, m, N, tol);

n_columns = size(H,2);

figure();

for j = 1:n_columns
 
    plot(real(all_eigenvalues), imag(all_eigenvalues), 'b*',...
    real(eigenvalues_krylov), imag(eigenvalues_krylov), 'g*',...
    real(H(:, j)), imag(H(:, j)), 'ro');
    legend('all eigenvalues','k desired','k found');
    drawnow
    % it seems like it freeze, but its just that the last eigenvalues are
    % all very close, so it seems like they are always the same
end

