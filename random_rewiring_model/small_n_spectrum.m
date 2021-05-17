function [ Evec_clean, Eval_clean, Evec, Eval, Eval_noise, spec_gap_clean, spec_gap_p] = small_n_spectrum(list, angle, list_f, angle_f, k_list, n, p)
%This is for testing how the spectral gaps change at different frequencies.
%This code only applies to the random rewiring model without small angle perturbation
% Input:
%   list: clean nearest neighbor list
%   angle: clean alignment angles
%   list_f: perturbed nearest neighbor list
%   angle_f: perturbed angles
%   k_list: list of frequencies
%
% Output: 
%   Evec_clean: the eigenvectors of the clean H^{(k)} matrices.
%   Eval_clean: the eigenvalues of the clean H^{(k)} matrices.
%   Evec: the eigenvectors of the noisy $H^{(k)}$ matrices. 
%   Eval: the eigenvalues of the noisy $H^{(k)}$ matrices. 
%   Eval_noise: the eigenvalues of the residual matrix
%   spec_gap_clean:
%   spec_gap_p:
%
% Zhizhen Zhao 04/2021

lk = length(k_list);
Evec_clean = cell(lk, 1);
Eval_clean = cell(lk, 1);
Evec = cell(lk, 1);
Eval = cell(lk, 1);
Eval_noise = cell(lk, 1);

Evec_noise =cell(lk, 1);
spec_gap_clean = zeros(lk, 1);
spec_gap_p = zeros(lk, 1);

for i = 1: lk
    [ sorted_eigval, sorted_eigvec, gap, W_clean ] = eig_decomp( list, angle, k_list(i), n );
    Eval_clean{i} = sorted_eigval;
    Evec_clean{i} = sorted_eigvec;
    spec_gap_clean(i) = gap;
    [ sorted_eigval, sorted_eigvec, gap, W_p ] = eig_decomp( list_f, angle_f, k_list(i), n );
    Eval{i} = sorted_eigval;
    Evec{i} = sorted_eigvec;
    spec_gap_p(i) = gap;
    N = W_p - p*W_clean; %the residual matrix in Eq. (6.6)
    [ u, d ] = eig(full(N)); %full spectrum
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    Eval_noise{i} = sorted_eigval;
    Evec_noise{i} = u(:, id);
end

end

function [ sorted_eigval, sorted_eigvec, spec_gap, W ] = eig_decomp( list, angle, k, n )

    W = sparse(list(:, 1), list(:, 2), exp(sqrt(-1)*k*angle*pi/180), n, n);
    W = W + W';
    D = sum(abs(W), 2);
    mk = 2*k + 1;
    %W = bsxfun(@times, 1./sqrt(D), W);
    %W = bsxfun(@times, 1./sqrt(D).', W);
    [ u, d ] = eig(full(W)); %compute the full spectrum
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    spec_gap = sorted_eigval(mk) - sorted_eigval(mk+1);
    sorted_eigvec = u(:, id);
    
end