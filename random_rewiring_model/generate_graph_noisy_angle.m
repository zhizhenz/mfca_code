function [ list, angle, list_f, angle_f, corr, angle_perturb ] = generate_graph_noisy_angle(q, p, threshold, kappa)

%This code generates initial graph under random rewiring model with small angular
%perturbation
% Input: 
%   q: quaternions, size 4\times n with n quaternions 
%   p: the probability to keep the clean edge
%   threshold: cutoff in the clean graph connection
%   kappa: parameter in the von Mises distribution. Kappa > 10^5 means no
%   small angular perturbation
%
% Output:
%   list: list of clean nearest neighbors
%   angle: list of clean alignment angles
%   list_f: list of nearest neighbors after perturbation
%   angle_f: list of noisy angles 
%   corr: clean correlation of the viewing angles
%   angle_perturb: \varepsilon small angular perturbation
%
% Zhizhen Zhao 04/2021

n = size(q, 2);
rot = q_to_rot(q);
v = squeeze(rot(3, :, :));
clear rot;
corr = v'*v;
corr2 = corr - diag(diag(corr));
list = find(corr2>threshold);
angle_perturb = 0;

[ I, J ] = ind2sub([n, n], list);
list = [ I, J ];
tmp = list(:, 1) - list(:, 2);
list = list(tmp>0, :);

angle = zeros(length(list), 1);
for i = 1:length(list)
    angle(i) = q_to_inplanerot(list(i, 1), list(i, 2), q);
end

corr = corr(1:100, 1:100); %clean correlation matrix

if p < 1
    ran_num = rand(length(list), 1);
    id1 = find(ran_num <= p);
    list_clean = list(id1, :);
    angle_clean = angle(id1);
    if kappa > 10^5
        angle_noisy = angle_clean;
    else
        [angle_perturb] = vmrand(0, kappa, length(id1), 1)*180/pi;
        angle_noisy = angle_clean + angle_perturb;
    end
    
    id2 = find(ran_num >p);
    lid2 = length(id2);
    list_p = list(id2, :);
    angle_p = 360*rand(lid2, 1);
    
    ran_num2 = rand(lid2, 1);
    id3 = find(ran_num2>0.5);
    id4 = find(ran_num2<=0.5);
    list1 = list_p(id3, :);
    list2 = list_p(id4, :);
    tmp_ind1 = randi([1, n-1], length(id3), 1);
    tmp_ind2 = randi([1, n-1], length(id4), 1);
    tmp_ind1(list1(:, 1) - tmp_ind1 <=0) = tmp_ind1(list1(:, 1) - tmp_ind1 <=0) + 1;
    list1(:, 2) = tmp_ind1;
    tmp_ind2(list2(:, 2) - tmp_ind2 <=0) = tmp_ind2(list2(:, 2) - tmp_ind2 <=0) + 1;
    list2(:, 1) = tmp_ind2;
    list_p = [list1; list2];
    
    list_p = sort(list_p, 2, 'descend');
    [ list_p, ia, ~ ]  = unique(list_p, 'rows');
    angle_p = angle_p(ia);
    
    list_f = [list_clean; list_p];
    angle_f = [angle_noisy; angle_p];
    [list_f, ia, ~] = unique(list_f, 'rows');
    angle_f = angle_f(ia);
else
    list_f = list;
    angle_f = angle;
end

