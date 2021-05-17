%Main jobscript for the multi-frequency class averaging on the
%probabilistic model detailed in Section 6

%%%%For Figures 4--9 in the paper
n = 10000; %number of quaternions 
initstate;
q = qrand(n); %quaternions
k_max = 50; %maximum frequency
k_list = [1:k_max];
k_list_scatter = [1, 5, 10];
k_list_hist = [1, 3, 5, 7];
k_list_eval = [1, 3, 5];
h = 0.92; %the theshold for the clean graph spherical cap size. h = 1-\cos a, where a is the opening angle for the spherical cap
kappa = 2*10^5; %variance for the von Mises distribution on the small angule perturbation. kappa >10^5 is treated as no angular perturbation.
nn = 50; %number of nearest neighbors
p_list = [ 1, 0.2, 0.1, 0.08 ]; %the random rewiring model, p is the percentage of the true connections, vary p to be 1, 0.2, 0.1, 0.08
eval_tmp = zeros(length(p_list), 3);
prop_threshold = 0.95;
k_max_joint = 20;

for j = 1:length(p_list)
    
    p = p_list(j);
    
    [ list, angle, list_f, angle_f, corr, angle_perturb ] = generate_graph_noisy_angle(q, p, h, kappa);
    
    [ feature, Evec, Eval, corr_ca, d, d_joint_g, d_joint_a, spec_gap ] = mfca(q, list_f, angle_f, k_max, nn);
    
    [ feature_single, corr_ca_single_2, d2 ] = ca(q, list_f, angle_f, k_list, nn);
    
    filename = sprintf('data_p%d_kappa%d_h%d_n%d_nn%d.mat', p*100, kappa, h*100, n, nn);
    save(filename, 'p', 'kappa', 'h', 'n', 'nn', 'corr', 'Evec', 'Eval', 'corr_ca', 'd', 'd_joint_g', 'd_joint_a', 'spec_gap', 'corr_ca_single_2', 'd2');
    
    %%Figure 4 bar plot of eigenvalues
    for i = 1:length(k_list_eval)
        k = k_list_eval(i);
        eval_tmp(:, i) = Eval{k}(1:19);
    end
    min_val = min(1-eval_tmp(:));
    max_val = max(1-eval_tmp(:));
    
    for i = 1:length(k_list_eval)
        k = k_list_eval(i);
        filename = sprintf('mfca_eval_n%d_p%d_kappa%d_h%d_k%d.fig', n, p*100, kappa, h*100, k);
        filename2 = sprintf('mfca_eval_n%d_p%d_kappa%d_h%d_k%d.png', n, p*100, kappa, h*100, k);
        figure; bar(1-Eval{k}(1:19));
        xlabel('$n - k + 1$', 'interpreter', 'latex');
        ylabel('$ 1 - \lambda_n^{(k)}$', 'interpreter', 'latex');
        set(gca, 'Fontsize', 22);
        if p == 1
            set(gca, 'Ylim', [min_val - 0.005, max_val + 0.005]);
        else
            set(gca, 'Ylim', [min_val - 0.001, max_val + 0.001]);
        end
        set(gca, 'XTick', [0:5:20]);
        saveas(gcf, filename);
        saveas(gcf, filename2);
    end
    
    %%%Figure 5 scatter plots
    tmp_corr = zeros(100^2, 1);
    for i = 1:length(k_list_scatter)
        k = k_list_scatter(i);
        filename = sprintf('mfca_scatter_n%d_p%d_kappa%d_h%d_k%d.fig', n, p*100, kappa, h*100, k);
        filename2 = sprintf('mfca_scatter_n%d_p%d_kappa%d_h%d_k%d.png', n, p*100, kappa, h*100, k);
        val = (corr(:) + 1).^k /2^k;
        tmp_corr = tmp_corr + log(corr_ca{k}(:));
        figure; scatter(val(:), corr_ca{k}(:));
        xlabel('$(\langle \pi(x_i), \pi(x_j) \rangle + 1)^k/2^k$', 'interpreter', 'latex');
        ylabel('$A^{(k)}_{ij}$', 'interpreter', 'latex');
        set(gca, 'Fontsize', 24);
        set(gca, 'Xlim', [0, 1]);
        set(gca, 'Ylim', [0, 1]);
        saveas(gcf, filename);
        saveas(gcf, filename2);
    end
    
    %%% Figure 6 scatter plots combined
    filename = sprintf('mfca_scatter_joint_n%d_p%d_kappa%d_h%d.fig', n, p*100, kappa, h*100);
    filename2 = sprintf('mfca_scatter_joint_n%d_p%d_kappa%d_h%d.png', n, p*100, kappa, h*100);
    val = corr(:);
    figure; scatter(val(:), tmp_corr);
    xlabel('$\langle \pi(x_i), \pi(x_j) \rangle$', 'interpreter', 'latex');
    ylabel('$\log(A^{\mathrm{All}}_{ij})$', 'interpreter', 'latex');
    set(gca, 'Fontsize', 24);
    saveas(gcf, filename);
    saveas(gcf, filename2);
       
    %%% Figure 7 histograms
    figure;
    xbin = linspace(0, 180, 100);
    interval = xbin(2) - xbin(1);
    for i = 1:length(k_list_hist)
        [ N1, X1 ] = hist(acosd(d{i}(:)), xbin);
        plot(xbin, N1/sum(N1)/interval, 'linewidth', 2);
        hold on;
    end
    [ N2, ~ ] = hist(acosd(d_joint_g{k_max_joint}(:)), xbin);
    plot(xbin, N2/sum(N2)/interval, 'linewidth', 2);
    xlabel('acos($\langle \pi(x_i), \pi(x_j) \rangle$)', 'interpreter', 'latex');
    ylabel('Proportion', 'interpreter', 'Latex');
    set(gca, 'Fontsize', 24);
    legend({'$k = 1$', '$k = 3$', '$ k = 5 $', '$ k = 7$', '$A^{All}$'}, 'interpreter', 'latex', 'location', 'best');
    legend boxoff;
    filename = sprintf('mfca_hist_n%d_p%d_kappa%d_h%d.fig', n, p*100, kappa, h*100);
    filename2 = sprintf('mfca_hist_n%d_p%d_kappa%d_h%d.png', n, p*100, kappa, h*100);
    set(gca, 'Xlim', [0, 180]);
    set(gca, 'XTick', [0:30:180]);
    saveas(gcf, filename);
    saveas(gcf, filename2);
    
    %%% Figure 8
    stats_single = zeros(k_max, 1);
    stats_g = zeros(k_max, 1);
    stats_a = zeros(k_max, 1);
    stats2 = zeros(k_max, 1);
    
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
    set(gca,'Ylim', [0.4, 1]);
    xlabel('$k$', 'interpreter', 'latex');
    ylabel('Proportion', 'interpreter', 'Latex');
    legend boxoff;
    
    saveas(gcf, filename);
    saveas(gcf, filename2);
    
end
