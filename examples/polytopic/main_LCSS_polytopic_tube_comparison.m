clear;

%% conservatism comparison with varying eps_A
% coverage = size of feasible domain/size of maximum robust invariant set

eps_A_list = 0.05:0.05:0.4;
eps_B = 0.1;
sigma_w = 0.1;

horizon = 10;

num_iter = length(eps_A_list);
poly_sls_mpc_data = cell(1, num_iter);
tube_mpc_data = cell(1, num_iter);

for ii = 1:num_iter
    eps_A = eps_A_list(ii);
    [mpc, Z_tube] = problem_generator(eps_A, eps_B, sigma_w, horizon);
    
    % the terminal set is chosen as the maximum robust invariant set
    terminal_set = mpc.terminal_constr;
    x0_set =  terminal_set.grid(15); 
    
    % coverage of polytopic SLS MPC
    [poly_sls_mpc_diags] = mpc_coverage(mpc, x0_set, 'polytopic_SLS_MPC');
    poly_sls_mpc_data{ii} = poly_sls_mpc_diags;
    
    % coverage of tube MPC
    [tube_mpc_diags] = mpc_coverage(mpc, x0_set, 'tube_MPC', Z_tube);
    tube_mpc_data{ii} = tube_mpc_diags;
end

poly_sls_mpc_coverage = zeros(1, num_iter);
tube_mpc_coverage = zeros(1, num_iter);
for ii = 1:num_iter
    poly_sls_mpc_coverage(ii) = poly_sls_mpc_data{ii}.coverage;
    tube_mpc_coverage(ii) = tube_mpc_data{ii}.coverage;
end

% load 'LCSS_conservatism_eps_A.mat' to generate figure in the paper.

figure; hold on;
plot(eps_A_list, poly_sls_mpc_coverage,'linewidth', 2, 'linestyle', '-');
plot(eps_A_list, tube_mpc_coverage, 'linewidth', 2, 'linestyle', '-.');
legend('poly-SLS-MPC', 'tube-MPC', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('$\epsilon_A$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('coverage', 'FontSize', 18, 'Interpreter', 'latex');
grid on

% save('LCSS_conservatism_eps_A.mat');

%% conservatism comparison with varying sigma_w
eps_A = 0.1;
eps_B = 0.1;
sigma_w_list = 0.1:0.1:0.7;

horizon = 10;

num_iter = length(sigma_w_list);
poly_sls_mpc_data_w = cell(1, num_iter);
tube_mpc_data_w = cell(1, num_iter);

for ii = 1:num_iter
    sigma_w = sigma_w_list(ii);
    [mpc, Z_tube] = problem_generator(eps_A, eps_B, sigma_w, horizon);
    
    % the terminal set is chosen as the maximum robust invariant set
    terminal_set = mpc.terminal_constr;
    x0_set =  terminal_set.grid(15); 
    
    % coverage of polytopic SLS MPC
    [poly_sls_mpc_diags] = mpc_coverage(mpc, x0_set, 'polytopic_SLS_MPC');
    poly_sls_mpc_data_w{ii} = poly_sls_mpc_diags;
    
    % coverage of tube MPC
    [tube_mpc_diags] = mpc_coverage(mpc, x0_set, 'tube_MPC', Z_tube);
    tube_mpc_data_w{ii} = tube_mpc_diags;
end

poly_sls_mpc_coverage = zeros(1, num_iter);
tube_mpc_coverage = zeros(1, num_iter);
for ii = 1:num_iter
    poly_sls_mpc_coverage(ii) = poly_sls_mpc_data_w{ii}.coverage;
    tube_mpc_coverage(ii) = tube_mpc_data_w{ii}.coverage;
end

% load 'LCSS_conservatism_sigma_w.mat' to generate figure in the paper

figure; hold on;
plot(sigma_w_list, poly_sls_mpc_coverage,'linewidth', 2, 'linestyle', '-');
plot(sigma_w_list, tube_mpc_coverage, 'linewidth', 2, 'linestyle', '-.');
legend('poly-SLS-MPC', 'tube-MPC', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('$\sigma_w$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('coverage', 'FontSize', 18, 'Interpreter', 'latex');
grid on

% save('LCSS_conservatism_sigma_w.mat');

%% calculate average solver time
% eps_A experiment
poly_sls_mpc_eps_A_runtime = 0;
poly_sls_mpc_eps_A_count = 0;

tube_mpc_eps_A_runtime = 0;
tube_mpc_eps_A_count = 0;
for ii = 1:8
   poly_sls_mpc_eps_A_runtime = poly_sls_mpc_eps_A_runtime  + poly_sls_mpc_data{ii}.running_time;
   poly_sls_mpc_eps_A_count = poly_sls_mpc_eps_A_count + size(poly_sls_mpc_data{ii}.x0_set, 1);
   tube_mpc_eps_A_runtime = tube_mpc_eps_A_runtime + tube_mpc_data{ii}.running_time;
   tube_mpc_eps_A_count = tube_mpc_eps_A_count + size(tube_mpc_data{ii}.x0_set, 1);
end

poly_sls_mpc_eps_A_avg_runtime = poly_sls_mpc_eps_A_runtime/poly_sls_mpc_eps_A_count;
tube_mpc_eps_A_avg_runtime = tube_mpc_eps_A_runtime/tube_mpc_eps_A_count;

% sigma_w experiment
poly_sls_mpc_w_runtime = 0;
poly_sls_mpc_w_count = 0;

tube_mpc_w_runtime = 0;
tube_mpc_w_count = 0;
for ii = 1:7
   poly_sls_mpc_w_runtime = poly_sls_mpc_w_runtime  + poly_sls_mpc_data_w{ii}.running_time;
   poly_sls_mpc_w_count = poly_sls_mpc_w_count + size(poly_sls_mpc_data_w{ii}.x0_set, 1);
   tube_mpc_w_runtime = tube_mpc_w_runtime + tube_mpc_data_w{ii}.running_time;
   tube_mpc_w_count = tube_mpc_w_count + size(tube_mpc_data_w{ii}.x0_set, 1);
end

poly_sls_mpc_w_avg_runtime = poly_sls_mpc_w_runtime/poly_sls_mpc_w_count;
tube_mpc_w_avg_runtime = tube_mpc_w_runtime/tube_mpc_w_count;

poly_avg_runtime = (poly_sls_mpc_eps_A_runtime + poly_sls_mpc_w_runtime)/(poly_sls_mpc_eps_A_count + poly_sls_mpc_w_count);
tube_avg_runtime = (tube_mpc_eps_A_runtime + tube_mpc_w_runtime)/(tube_mpc_eps_A_count + tube_mpc_w_count);



