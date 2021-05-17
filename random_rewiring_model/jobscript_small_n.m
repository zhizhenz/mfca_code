%Main jobscript for the multi-frequency class averaging on the
%probabilistic model detailed in Section 6

%%%%For Figures 2--3 in the paper
n = 1000; %number of quaternions
%initstate;
q = qrand(n); %quaternions
k_max = 20; %maximum frequency
%k_list = [1, 4, 8, 20];
k_list = 1:12;
h = 0.8; %the theshold for the clean graph spherical cap size. h = 1-\cos a, where a is the opening angle for the spherical cap
kappa = 2*10^5; %variance for the von Mises distribution on the small angule perturbation. kappa >10^5 is treated as no angular perturbation.
nn = 50; %number of nearest neighbors
p_list = [ 0.5, 0.3, 0.15 ]; %the random rewiring model, p is the percentage of the true connections, vary p to be 0.5, 0.3 and 0.15
prop_threshold = 0.9;

for j = 1:length(p_list)
    
    p = p_list(j);
    [ list, angle, list_f, angle_f, corr, angle_perturb ] = generate_graph_noisy_angle(q, p, h, kappa);
    
    [ Evec_clean, Eval_clean, Evec, Eval, Eval_noise, spec_gap_clean, spec_gap_p] = small_n_spectrum(list, angle, list_f, angle_f, k_list, n, p);
    
    filename = sprintf('smalln_data_eigval_p%d_kappa%d_h%d_n%d_nn%d.mat', p*100, kappa, h*100, n, nn);
    save(filename, 'Eval_clean', 'Eval', 'Eval_noise', 'spec_gap_clean', 'spec_gap_p');
    
    xbin = linspace(-20, 60, 60);
    
    for i = 1:length(k_list)
        
        [ N1, X1 ] = hist(Eval{i}(:), xbin);
        [ N2, X2 ] = hist(Eval_noise{i}(:), xbin);
        
        figure; bar(X1, N1);
        filename1 = sprintf('spec_n1000_unormalized_p%d_k%d.fig', p*100, k_list(i));
        filename2 = sprintf('spec_n1000_unormalized_p%d_k%d.png', p*100, k_list(i));
        set(gca, 'Fontsize', 20);
        set(gca, 'Xlim', [-20, 60]);
        set(gca, 'Ylim', [0, 55]);
        xlabel('$\lambda$', 'interpreter', 'latex');
        ylabel('Counts', 'interpreter', 'latex');
        saveas(gcf, filename1);
        saveas(gcf, filename2);
        
        figure; bar(X2, N2);
        filename1 = sprintf('spec_n1000_noise_p%d_k%d.fig', p*100, k_list(i));
        filename2 = sprintf('spec_n1000_noise_p%d_k%d.png', p*100, k_list(i));
        set(gca, 'Fontsize', 20);
        set(gca, 'Xlim', [-20, 60]);
        set(gca, 'Ylim', [0, 55]);
        xlabel('$\lambda$', 'interpreter', 'latex');
        ylabel('Counts', 'interpreter', 'latex');
        saveas(gcf, filename1);
        saveas(gcf, filename2);
        
    end   
       
    
    [ feature, Evec, Eval, corr_ca, d, d_joint_g, d_joint_a, spec_gap ] = mfca(q, list_f, angle_f, k_max, nn);
    
    %%% For Figure 3
    stats = zeros(k_max, 1);

    
    for i = 1:k_max
        stats(i) = length(find(d{i}(:)>prop_threshold))/(nn*n);
    end
    
    k_list2 = 1:k_max;
    filename = sprintf('CA_performance_n_%d_vark_p%d_th%d_kappa%d_nn%d_h%d.fig',n, p*100, prop_threshold*100, kappa, nn, h*100);
    filename2 = sprintf('CA_performance_n_%d_vark_p%d_th%d_kappa%d_nn%d_h%d.png', n, p*100, prop_threshold*100, kappa, nn, h*100);
    figure; plot(k_list2, stats, 'linewidth', 2);
    set(gca, 'Fontsize', 20);
    xlabel('$k$', 'interpreter', 'latex');
    ylabel('Proportion', 'interpreter', 'Latex');
    saveas(gcf, filename);
    saveas(gcf, filename2);
    
    filename = sprintf('smalln_data_mfca_p%d_kappa%d_h%d_n%d_nn%d_th%d.mat', p*100, kappa, h*100, n, nn, prop_threshold*100);
    save(filename, 'Evec', 'Eval', 'd', 'd_joint_g', 'd_joint_a', 'stats', 'prop_threshold');

    
end