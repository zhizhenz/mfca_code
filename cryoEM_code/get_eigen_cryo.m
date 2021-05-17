%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build initial nearest neighbor graph and eigen-decompositon of $H^{(k)}$
% Inputs:
%   angle: rotation alignment between nearest neighbors following the list
%   list: list of nearest neighbors 
%   eigen_num: m_k for k = 1,...,k_max
% Outputs:
%   Eval: top $m_k$ eigenvalues of $H^{(k)}$ for k = 1,...,k_max
%   Evec: corresponding eigenvectors
%
% Yifeng Fan, 2021/04/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Eval, Evec ] = get_eigen_cryo( angle, list, refl, eigen_num )

%%% Build $H^{(k)}$ %%%
n = size(list, 1);
knn = size(list, 2);
list = [ repmat([1:n]', knn, 1), list(:)];
refl = refl(:);
angle = angle(:);

list_nrefl = list(refl == 1,:);
list_refl = list(refl == 2,:);
angle_nrefl = angle(refl == 1,:);
angle_refl = angle(refl == 2,:);

[ list_nrefl, id ] = sort(list_nrefl,2,'ascend'); 
angle_nrefl(id(:,1) == 2) = -angle_nrefl(id(:,1) == 2); 
[ list_nrefl, id2 ] = unique(list_nrefl,'rows'); 
angle_nrefl = angle_nrefl(id2);
list_nrefl = [ list_nrefl; list_nrefl+n ]; 
angle_nrefl = [-angle_nrefl; angle_nrefl ];

[ list_refl, ~ ] = sort(list_refl,2,'ascend'); 
[ list_refl, id2 ] = unique(list_refl,'rows'); 
angle_refl = angle_refl(id2);
list_refl = [ list_refl(:,1), list_refl(:,2)+n; list_refl(:,2), list_refl(:,1)+n ];
angle_refl = [180-angle_refl;180-angle_refl ];

I_1tok = [list_nrefl(:,1); list_refl(:,1)];
J_1tok = [list_nrefl(:,2); list_refl(:,2)];
angle_1tok = [ angle_nrefl; angle_refl ];

%%% Eigen-decomposition of $H^{(k)}$ for each k = 1,...,k_max %%%
parfor i = 1: numel(eigen_num)
    
    % Construct $H^{(k)}$
    H_k = sparse(I_1tok, J_1tok, exp(sqrt(-1)*(i)*angle_1tok*pi/180), 2*n, 2*n);
    H_k = H_k + H_k';
    H_k = abs(H_k);
    D = sum(H_k, 2);
    H_k = bsxfun(@times, 1./sqrt(D), H_k);
    H_k = bsxfun(@times, 1./sqrt(D).', H_k);
    
    % Eigen-decomposition
    [ u, d ] = eigs(H_k, eigen_num(i));
    [ sorted_eigval, id ] = sort(real(diag(d)), 'descend');
    sorted_eigvec = u(:, id);
    sorted_eigval(isnan(sorted_eigval)) = 0;
    sorted_eigvec(:,isnan(sorted_eigval)) = 0; 
    Evec{i} = sorted_eigvec;
    Eval{i} = sorted_eigval;
    
end


end

