%Main jobscript for the multi-frequency class averaging on the
%probabilistic model detailed in Section 6

%%%%For Figures 10 in the paper
n = 10000; %number of quaternions 
initstate;
q = qrand(n); %quaternions
k_max = 10; %maximum frequency
k_list = [1:k_max];
k_list_scatter = [1, 5, 10];
k_list_hist = [1, 3, 5, 7];
l = 2*k_max + 1;
h = 0.7; %the theshold for the clean graph spherical cap size. h = 1-\cos a, where a is the opening angle for the spherical cap
%kappa = 2*10^5; %variance for the von Mises distribution on the small angule perturbation. kappa >10^5 is treated as no angular perturbation.
nn = 50; %number of nearest neighbors
p = 0.08; %the random rewiring model, p is the percentage of the true connections, vary p to be 1, 0.2, 0.1, 0.08
kappa_list = [ 2*10^5, 500, 64];

for j = 1:length(kappa_list)
    
    kappa = kappa_list(j);
    
    [ list, angle, list_f, angle_f, corr, angle_perturb ] = generate_graph_noisy_angle(q, p, h, kappa);
    
    [ feature, Evec, Eval, corr_ca, d, d_joint_g, d_joint_a, spec_gap ] = mfca(q, list_f, angle_f, k_max, nn);
    
    [ feature_single, corr_ca_single_2, d2 ] = ca(q, list_f, angle_f, k_list, nn);
    
    filename = sprintf('data_p%d_kappa%d_h%d_n%d_nn%d.mat', p*100, kappa, h*100, n, nn);
    save(filename, 'p', 'kappa', 'h', 'n', 'nn', 'corr', 'Evec', 'Eval', 'corr_ca', 'd', 'd_joint_g', 'd_joint_a', 'spec_gap', 'corr_ca_single_2', 'd2');
    
    %%% Figure 10
    stats_single = zeros(k_max, 1);
    stats_g = zeros(k_max, 1);
    stats_a = zeros(k_max, 1);
    stats2 = zeros(k_max, 1);
    prop_threshold = 0.95;
    
    for i = 1:k_max
        stats_single(i) = length(find(d{i}(:)>prop_threshold))/(nn*n);
        stats_g(i) = length(find(d_joint_g{i}(:)>prop_threshold))/(nn*n);
        stats_a(i) = length(find(d_joint_a{i}(:)>prop_threshold))/(nn*n);
        stats2(i) = length(find(d2{i}(:)>prop_threshold))/(nn*n);
    end
    
    filename = sprintf('CA_performance_n_%d_vark_p%d_th%d_kappa%d_nn%d_h%d.fig',n, p*100, prop_threshold*100, kappa, nn, h*100);
    filename2 = sprintf('CA_performance_n_%d_vark_p%d_th%d_kappa%d_nn%d_h%d.png', n, p*100, prop_threshold*100, kappa, nn, h*100);
    figure; plot(k_list, stats_single, k_list, stats_g, k_list, stats_a, k_list, stats2, 'linewidth', 2);
    legend({'$A^{(k)}$', '$A^{All}$', '$S^{All}$', '$B^{(k)}$'}, 'interpreter', 'latex', 'location', 'best');
    set(gca, 'Fontsize', 20);
    set(gca, 'Xlim', [1, k_max]);
    xlabel('$k$', 'interpreter', 'latex');
    ylabel('Proportion', 'interpreter', 'Latex');
    legend boxoff;
    
    saveas(gcf, filename);
    saveas(gcf, filename2);
    
end
