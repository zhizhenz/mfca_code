%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for reproducing the cryo-EM experimental result in 'Representation 
% Theoretic Patterns in Multi-Frequency Class Averaging for Three-Dimensional 
% Cryo-Electron Microscopy'
% 
% Please download ASPIRE package before running the code.
%
% Yifeng Fan, 2021/04/19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Please run the following code at first %%%
clear
close all
addpath /home/aspire;       %add the ASPIRE package to the path
initpath
rng('default');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters setting %%%
SNR = 0.008;                % SNR of the noisy image
knn_ini = 50;               % Number of nearest neighbors  
eigen_num = 2*[1:20]+1;     % m_k for k = 1,...k_max


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preprocessing $$$
load /home/clean.mat;          % Clean projection images
load /home/quaternion.mat;     % Quaternions of the projection images
q = q(:,1:10000); 

images = struct();
[images.noisy, noise_v_r] = addnoise_v6(real(projs), SNR);      % Add while Gaussian noise

% Support estimate
l = size(images.noisy, 1);
n = size(images.noisy, 3);
R = floor(l/2);
[ x, y ] = meshgrid(-R:R, -R:R);
r = sqrt(x.^2 + y.^2);
threshold = 0.999;
[ c, R ] = support_estimate(images.noisy, threshold);
support.c = c;
support.R = R;
support.R = 60;

% Multiple nodes
delete(gcp)
parpool('local',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial nearest neighbor search $$$
log_message('Start initial nearest neighbor search')
[ sPCA_data, ~, ~, ~, ~ ] =  data_sPCA_update(images.noisy, noise_v_r, support);
[ class.sPCA, refl.sPCA, rot.sPCA, ~,  ~ ] = Initial_classification_FD_update(sPCA_data, knn_ini, 0);
log_message('Finish initial nearest neighbor search')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eigen-decomposition of $H^{(k)}$ %%%
log_message('Start eigen-decomposition');
[ Eval, Evec ] = get_eigen_cryo(rot.sPCA, class.sPCA, refl.sPCA, eigen_num);
log_message('Finish eigen-decomposition');

% % Plot of the spectrum in Fig.12
% for i = [1, 3, 5]
%     figure; bar(Eval{i})
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MFCA, Compute Affinity measure $A^{all}$%%%

aff_mfca = mfca_k(Evec);    

%%% Scatter plot in Fig.13 and Fig.14 %%%

% Using $A^{(all)}$
aff_tmp = log(aff_mfca(1:100,1:100));
rot_tmp = q_to_rot(q);
v = squeeze(rot_tmp(3, :, :)); 
tmp_x = v(:,1:100)'*v(:,1:100);
figure; scatter(tmp_x(:), aff_tmp(:), 'o')

% Using single frequency $A^{(k)}$
freq_range = [1,5,10]; 
tmp_rot = q_to_rot(q);
v = squeeze(tmp_rot(3, :, :)); 
for num = 1:numel(freq_range)
    i = freq_range(num);
    aff_single = mfca_k(Evec(i));
    aff_single = aff_single(1:100,1:100);
    tmp_x = ((v(:,1:100)'*v(:,1:100) + 1).^i)./(2^i);
    figure; scatter(tmp_x(:), abs(aff_single(:)), 'o')
end

%%% Histogram in Fig.15 %%%

% Using single frequency
f_range = [1,3,5];
legend_text = cell(1, numel(f_range)+1);
figure;
for num = 1:numel(f_range)
    i = f_range(num);
    tmp_aff = mfca_k(Evec(i));
    tmp_aff = tmp_aff - diag(diag(tmp_aff));
    knn = 50;
    [~, id_tmp] = sort(tmp_aff(1:n,:), 2, 'descend');
    class_tmp = id_tmp(:,1:knn);
    refl_tmp = (class_tmp > n) + 1;
    class_tmp(class_tmp > n) = class_tmp(class_tmp > n) - n;
    [ e_c_tmp, ~ ] = check_simulation_results(class_tmp, refl_tmp, zeros(size(class_tmp)), q);
    [ x_c_tmp, y_c_tmp] = hist(acos(e_c_tmp(:))*180/pi,0:180);
    plot(y_c_tmp,x_c_tmp/sum(x_c_tmp))
    hold on
    legend_text{1,num} = ['k = ', num2str(i)];
end

% Using $A^{(all)}$
tmp_aff = aff_mfca - diag(diag(aff_mfca));
knn = 50;
[~, id_tmp] = sort(tmp_aff(1:n,:), 2, 'descend');
class_tmp = id_tmp(:,1:knn);
refl_tmp = (class_tmp > n) + 1;
class_tmp(class_tmp > n) = class_tmp(class_tmp > n) - n;
[ e_c_tmp, ~ ] = check_simulation_results(class_tmp, refl_tmp, zeros(size(class_tmp)), q);
[ x_c_tmp, y_c_tmp] = hist(acos(e_c_tmp(:))*180/pi,0:180);
plot(y_c_tmp,x_c_tmp/sum(x_c_tmp))
legend_text{1, numel(f_range)+1} = 'MFCA';
legend(legend_text);

