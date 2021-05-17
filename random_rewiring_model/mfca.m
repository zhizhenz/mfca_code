function [ feature, Evec, Eval, corr_ca, dist, dist_joint_g, dist_joint_a, spec_gap ] = mfca(q, list, angle, k_max, nn)

%This code for mfca
% Input: 
%   q: quaternions, size 4\times n with n quaternions 
%   list: nearest neighbor list
%   angle: pairwise alignment angle
%   k_max: maximum angular frequency cutoff
%   nn: number of nearest neighbors
%
% Output:
%   feature: list of clean nearest neighbors
%   Evec: the eigenvectors of \widetilde{H^{(k)} at different k
%   Eval: the eigenvalues of \widetilde{H^{(k)} at different k
%   corr_ca: the 100x100 correlation matrix using the mfca features
%   dist: angular distance between identified nearest neighbors using
%   A^{(k)}
%   dist_joint_g: angular distance between identified nearest neighbors using
%   A^{(All)}
%   dist_joint_a: angular distance between identified nearest neighbors using
%   S^{(All)}
%   spec_gap: \varepsilon small angular perturbation
%
% Zhizhen Zhao 04/2021

n = max(list(:));
Evec = cell(k_max, 1);
Eval = cell(k_max, 1);
feature = cell(k_max, 1);
corr_ca = cell(k_max, 1);
dist = cell(k_max, 1);
dist_joint_g = cell(k_max, 1);
dist_joint_a = cell(k_max, 1);

tmp_joint_g = ones(n);
tmp_joint_a = zeros(n);
spec_gap = zeros(k_max, 1);

for i = 1: k_max
    W = sparse(list(:, 1), list(:, 2), exp(sqrt(-1)*i*angle*pi/180), n, n);
    W = W + W';
    D = sum(abs(W), 2);
    mk = 2*i + 1;
    W = bsxfun(@times, 1./sqrt(D), W);
    W = bsxfun(@times, 1./sqrt(D).', W);
    [ u, d ] = eigs(W, max(30, 2*mk)); %top eigenvectors
    %[ u, d ] = eigs(W, 50); %top 2k + 1 eigenvectors
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    spec_gap(i) = sorted_eigval(mk) - sorted_eigval(mk+1);
    sorted_eigvec = u(:, id);
    norm_row = sqrt(sum(abs(sorted_eigvec(:, 1:mk)).^2, 2));
    feature{i} = bsxfun(@times, sorted_eigvec(:, 1:mk), 1./norm_row);
    Eval{i} = sorted_eigval;
    Evec{i} = sorted_eigvec;
    
    if n < 2000
        tmp_corr = zeros(n, n, i);
        tmp_corr(:, :, i) = feature{i}*feature{i}';
        tmp = abs(tmp_corr(:, :, i));
    else
        tmp = abs(feature{i}*feature{i}');
    end
    
    tmp_joint_g = tmp_joint_g.*tmp;
    
    tmp_joint_a = tmp_joint_a + 2*(tmp.^(1/i)) -1;
    
    corr_ca{i} = tmp(1:100, 1:100);
    
    tmp = tmp - diag(ones(n, 1));
    [ ~, id ] = sort(tmp, 2, 'descend');
    
    tmp_joint_g = tmp_joint_g - diag(ones(n, 1));
    [ ~, id_joint_g ] = sort(tmp_joint_g, 2, 'descend');
    
    tmp_joint_a = tmp_joint_a - diag(diag(tmp_joint_a));
    [ ~, id_joint_a ] = sort(tmp_joint_a, 2, 'descend');

    
    dist{i} = zeros(n, nn);
    dist_joint_g{i} = zeros(n, nn);
    dist_joint_a{i} = zeros(n, nn);
    
    for s = 1:n
        for t = 1:nn
            dist{i}(s, t) = q_to_d(s, id(s, t), q);
            dist_joint_g{i}(s, t) = q_to_d(s, id_joint_g(s, t), q);
            dist_joint_a{i}(s, t) = q_to_d(s, id_joint_a(s, t), q);
        end
    end
end

end
