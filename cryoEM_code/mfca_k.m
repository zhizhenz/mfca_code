%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the MFCA affinity measurement $A^{(all)}$
% Inputs:
%   Evec: eigenvectors of $H^{(k)}$ for k = 1,...,k_max
% Outputs:
%   affinity: MFCA affinity $A^{(all)}$
%
% Yifeng Fan, 2021/04/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ affinity ] = mfca_k( Evec )

k_max = size(Evec,2);
n = size(Evec{1},1);
affinity = ones(n,n);
for i = 1:k_max
    evec = Evec{i};
    evec_norm = sqrt(sum(abs(evec).^2,2));
    Evec{i} = diag(1./evec_norm)*evec;
    affinity = affinity.*abs(Evec{i}*Evec{i}');
end

end

