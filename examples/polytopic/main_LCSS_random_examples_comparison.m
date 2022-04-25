% compare the feasible domain of polytopic SLS MPC and tube MPC over
% randomly generated systems
clear;

num_expr = 100;
poly_sls_mpc_data = cell(1, num_expr);
tube_mpc_data = cell(1, num_expr);

eps_A = 0.2; eps_B = 0.1; sigma_w = 0.1; horizon = 10;

for ii = 1:num_expr
    [mpc, Z_tube] = random_problem_generator(eps_A, eps_B, sigma_w, horizon);
    
    % grid search of initial conditions over the state space
    x0_set =  mpc.state_constr.grid(15); 
    
    % coverage of polytopic SLS MPC
    [poly_sls_mpc_diags] = mpc_coverage(mpc, x0_set, 'polytopic_SLS_MPC');
    poly_sls_mpc_data{ii} = poly_sls_mpc_diags;
    
    % coverage of tube MPC
    [tube_mpc_diags] = mpc_coverage(mpc, x0_set, 'tube_MPC', Z_tube);
    tube_mpc_data{ii} = tube_mpc_diags;
    
%     save('conservatism_comparison_random_example.mat');
end


%% plot the coverage of polytopic SLS MPC in an ascending order 
% load 'conservatism_comparison_random_example.mat' to regenerate the
% figure in the paper. 

poly_sls_mpc_coverage = zeros(1, num_expr);
tube_mpc_coverage = zeros(1, num_expr);
for ii = 1:num_expr
    poly_sls_mpc_coverage(ii) = poly_sls_mpc_data{ii}.coverage;
    tube_mpc_coverage(ii) = tube_mpc_data{ii}.coverage;
end

[~, ind] = sort(poly_sls_mpc_coverage);

figure; hold on;
plot(poly_sls_mpc_coverage(ind),'linewidth', 2, 'linestyle', '-');
plot(tube_mpc_coverage(ind), 'linewidth', 2, 'linestyle', '-.');
legend('poly-SLS-MPC', 'tube-MPC', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('example number', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('coverage', 'FontSize', 18, 'Interpreter', 'latex');
grid on

% calculate total runtime
poly_sls_mpc_runtime = 0;
tube_mpc_runtime = 0;

for ii = 1:num_expr
    poly_sls_mpc_runtime = poly_sls_mpc_runtime + poly_sls_mpc_data{ii}.running_time;
    tube_mpc_runtime = tube_mpc_runtime + tube_mpc_data{ii}.running_time;
end
