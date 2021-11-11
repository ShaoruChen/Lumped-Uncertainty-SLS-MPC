clear;

data = load('mpc_data.mat');
mpc = data.mpc;
Z_inv = data.Z_inv;

horizons = 1:15;
mpc.x0 = [1;0];
solver_time_seq = [];
for ii =progress(1:length(horizons))
    horizon = horizons(ii);
    mpc.horizon = horizon;
    
[aug_df_sol] = mpc.SolveAugDistFeedbackSLSMPC('value');
if aug_df_sol.status ~= 0 
    aug_df_time = nan;
else
    aug_df_time = aug_df_sol.solver_time;
end

[unif_df_sol] = mpc.SolveUniformDistFeedbackMPC('value', 2);
if unif_df_sol.status ~= 0
    unif_df_time = nan;
else
    unif_df_time = unif_df_sol.solver_time;
end

[tube_sol] = mpc.SolveTubeMPC(Z_inv, 'value', 2);
if tube_sol.diagnostics.problem~= 0
    tube_time = nan;
else
    tube_time = tube_sol.solver_time;
end

[grid_sls_sol] = mpc.SolveSLSMPCAuto();
if grid_sls_sol.status ~= 0
    grid_sls_time = nan;
else
    grid_sls_time = grid_sls_sol.solver_time;
end
solver_time = [aug_df_time unif_df_time tube_time grid_sls_time]';
solver_time_seq = [solver_time_seq solver_time];
end

% save('data\solver_time_horizon_small.mat');


%% plot
aug_sls_time = solver_time_seq(1,:);
unif_df_time = solver_time_seq(2,:);
tube_time = solver_time_seq(3, :);
naive_sls_time = solver_time_seq(4,:);

figure; hold on; grid on;
plot(horizons, aug_sls_time, '-', 'LineWidth', 2);
plot(horizons, unif_df_time, '-.', 'LineWidth', 2);
plot(horizons, tube_time, ':', 'LineWidth', 2);
plot(horizons, naive_sls_time, '--', 'LineWidth', 2	);

set(gca,'yscale','log')
legend({'aug-SLS', 'unif-df', 'tube', 'grid-SLS'}, 'Interpreter', 'latex', 'FontSize', 16);
xlabel('horizon $T$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('solver time [sec]', 'Interpreter', 'latex', 'FontSize', 18);

