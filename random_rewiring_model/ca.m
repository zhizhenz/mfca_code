function [ feature, corr_ca, d1] = ca(q, list, angle, klist, nn)

%This code for single frequency class averaging
% Input: 
%   q: quaternions, size 4\times n with n quaternions 
%   list: nearest neighbor list
%   angle: pairwise alignment angle
%   klist: the list of angular frequencies
%   nn: number of nearest neighbors
%
% Output:
%   feature: list of clean nearest neighbors
%   corr_ca: the 100x100 correlation matrix using the ca features
%   d: angular distance between identified nearest neighbors using
%   B^{(k)}
%
% Zhizhen Zhao 04/2021

n = max(list(:));
mk = 2*max(klist) + 1;

corr_ca = cell(length(klist), 1);
feature = cell(length(klist), 1);
d1 = cell(length(klist), 1);

W = sparse(list(:, 1), list(:, 2), exp(sqrt(-1)*angle*pi/180), n, n);
W = W + W';
%W = W + diag(diag(ones(n, 1)));
D = sum(abs(W), 2);
W = bsxfun(@times, 1./sqrt(D), W);
W = bsxfun(@times, 1./sqrt(D).', W);
[ u, d ] = eigs(W, mk + 10); %top mk eigenvectors
[ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
sorted_eigvec = u(:, id);

for i = 1:length(klist)
    norm_row = sqrt(sum(abs(sorted_eigvec(:, 1:(2*(klist(i))+1))).^2, 2));
    feature{i} = bsxfun(@times, sorted_eigvec(:, 1:(2*(klist(i))+1)), 1./norm_row);
    tmp = abs(feature{i}*feature{i}');
    corr_ca{i} = tmp(1:100, 1:100);
    tmp = tmp - diag(ones(n, 1));
    [ sorted_tmp, id ] = sort(tmp, 2, 'descend');
    d1{i} = zeros(n, nn);
    for s = 1:n
        for t = 1:nn
            d1{i}(s, t) = q_to_d(s, id(s, t), q);
        end
    end
end

end
