clear;

data = load('mpc_data.mat');
mpc = data.mpc;
Z= data.Z_inv;

mpc.terminal_constr = [];

X_grid = mpc.state_constr;
x0_set = X_grid.grid(15);

%% test 1: varying epsilon_A
eps_A_list = 0.05:0.01:0.25;
eps_B = 0.1;
sigma_w = 0.1;

methods = {'aug_sls', 'tube', 'unif_df', 'grid_sls'};
mpc.eps_B = eps_B; mpc.sigma_w = sigma_w;

result = struct;
for kk = 3:4
    method = methods{kk};
    fprintf('Method: %s', method);
    diags_cell = cell(1, length(eps_A_list));
    for ii = 1:length(eps_A_list)
        eps_A = eps_A_list(ii);
        mpc.eps_A = eps_A;
        [diags] = MPC_feasibility_check(mpc, x0_set, method, Z);
        if isempty(diags.feasible_set)
            break;
        else
            diags_cell{ii} = diags;
        end
    end
    result.(method) = diags_cell;
    save('data\feasibility_varying_params_temp_0726.mat');
end

%% plot the coverage
aug_sls_count = zeros(1, length(eps_A_list));
for ii = 1:length(eps_A_list)
   if ~isempty(result.aug_sls{ii}) && ~isempty(result.aug_sls{ii}.feasible_set)
        aug_sls_count(ii) = size(result.aug_sls{ii}.feasible_set, 1);
   else
       break;
   end
end

unif_df_count = zeros(1, length(eps_A_list));
for ii = 1:length(eps_A_list)
   if ~isempty(result.unif_df{ii}) && ~isempty(result.unif_df{ii}.feasible_set)
        unif_df_count(ii) = size(result.unif_df{ii}.feasible_set, 1);
   else
       break;
   end
end


tube_count = zeros(1, length(eps_A_list));
for ii = 1:length(eps_A_list)
   if ~isempty(result.tube{ii}) && ~isempty(result.tube{ii}.feasible_set)
        tube_count(ii) = size(result.tube{ii}.feasible_set, 1);
   else
       break;
   end
end

grid_sls_count = zeros(1, length(eps_A_list));
for ii = 1:length(eps_A_list)
   if ~isempty(result.grid_sls{ii}) && ~isempty(result.grid_sls{ii}.feasible_set)
        grid_sls_count(ii) = size(result.grid_sls{ii}.feasible_set, 1);
   else
       break;
   end
end

total_num = size(x0_set, 1);
figure; hold on
plot(eps_A_list,aug_sls_count/total_num,'linewidth', 2, 'linestyle', '-');
plot(eps_A_list,unif_df_count/total_num, 'linewidth', 2, 'linestyle', '-.');
plot(eps_A_list,tube_count/total_num, 'linewidth', 2, 'linestyle', '--');
plot(eps_A_list,grid_sls_count/total_num, 'linewidth', 2, 'linestyle', ':');
legend('aug-SLS', 'unif-df', 'tube', 'grid-SLS', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('$\epsilon_A$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('coverage', 'FontSize', 18, 'Interpreter', 'latex');
grid on

%%

% mpc.eps_A = 0.1; mpc.eps_B = 0.1; mpc.sigma_w = 0.3;
% [diags] = MPC_feasibility_check(mpc, x0_set, 'tube', Z);
% 
% % [tube_sol] = mpc.SolveTubeMPC(Z, 'value',2);
% % [traj_set] = mpc.SimulateTrajwithTubeMPC( tube_sol, 10);
% % PlotFcns.plot_tube_MPC_traj(mpc, tube_sol, traj_set, 'tube mpc');
% 
% figure;
% scatter(x0_set(:,1), x0_set(:,2));
% hold on
% if ~isempty(diags.feasible_set)
%     scatter(diags.feasible_set(:,1), diags.feasible_set(:,2),'r*');
% end



