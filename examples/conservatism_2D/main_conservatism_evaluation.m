clear;

data = load('mpc_data.mat');
mpc = data.mpc;
Z_inv = data.Z_inv;

X_grid = mpc.terminal_constr;
x0_set = X_grid.grid(20);

%% uniform disturbance feedback MPC evaluation
num_points = size(x0_set, 1);
feasibleSet =[]; infeasibleSet = [];
start_time = tic;
verbose = 0;
[unif_df_mpc_optimizer] = mpc.SolveUniformDistFeedbackMPC('optimizer', verbose);

start_time = tic;
for ii = progress(1:num_points)
    init_x = x0_set(ii, :)';
    [sol_value, errorcode] = unif_df_mpc_optimizer(init_x);
    
    if errorcode == 0
        feasibleSet = [feasibleSet; init_x'];
    else
        infeasibleSet = [infeasibleSet; init_x'];
    end
end
uniform_df_MPC_running_time = toc(start_time);

unif_df_mpc_diags = struct;
unif_df_mpc_diags.feasible_set = feasibleSet;
unif_df_mpc_diags.infeasible_set = infeasibleSet;
unif_df_mpc_diags.running_time = uniform_df_MPC_running_time;
unif_df_mpc_diags.x0_set = x0_set;

figure;
feasibleSet = unif_df_mpc_diags.feasible_set;
infeasibleSet = unif_df_mpc_diags.infeasible_set;
scatter(x0_set(:,1), x0_set(:,2));
hold on
if ~isempty(feasibleSet)
    scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
end
title('uniform disturbance feedback');
xlabel('x_1'); ylabel('x_2'); 

%% augmented disturbance feedback MPC evaluation
num_points = size(x0_set, 1);
feasibleSet =[]; infeasibleSet = [];
start_time = tic;
opt = struct;
opt.solver = 'mosek'; opt.verbose = 0;
[aug_df_mpc_optimizer] = mpc.SolveAugDistFeedbackSLSMPC('optimizer', opt);

start_time = tic;
for ii = progress(1:num_points)
    init_x = x0_set(ii, :)';
    [sol_value, errorcode] = aug_df_mpc_optimizer(init_x);
    
    if errorcode == 0
        feasibleSet = [feasibleSet; init_x'];
    else
        infeasibleSet = [infeasibleSet; init_x'];
    end
end
aug_df_MPC_running_time = toc(start_time);

aug_df_mpc_diags = struct;
aug_df_mpc_diags.feasible_set = feasibleSet;
aug_df_mpc_diags.infeasible_set = infeasibleSet;
aug_df_mpc_diags.running_time = aug_df_MPC_running_time;
aug_df_mpc_diags.x0_set = x0_set;


figure;
feasibleSet = aug_df_mpc_diags.feasible_set;
infeasibleSet = aug_df_mpc_diags.infeasible_set;
scatter(x0_set(:,1), x0_set(:,2));
hold on
if ~isempty(feasibleSet)
    scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
end
title('augmented disturbance feedback');
xlabel('x_1'); ylabel('x_2'); 

%% tube MPC evaluation
num_points = size(x0_set, 1);
feasibleSet =[]; infeasibleSet = [];
start_time = tic;

verbose = 0;
Z_inv = data.Z_inv;
[tube_mpc_optimizer] = mpc.SolveTubeMPC(Z_inv, 'optimizer', verbose);

start_time = tic;
for ii = progress(1:num_points)
    init_x = x0_set(ii, :)';
    [sol_value, errorcode] = tube_mpc_optimizer(init_x);
    
    if errorcode == 0
        feasibleSet = [feasibleSet; init_x'];
    else
        infeasibleSet = [infeasibleSet; init_x'];
    end
end
tube_MPC_running_time = toc(start_time);


tube_mpc_diags = struct;
tube_mpc_diags.feasible_set = feasibleSet;
tube_mpc_diags.infeasible_set = infeasibleSet;
tube_mpc_diags.running_time = tube_MPC_running_time;
tube_mpc_diags.x0_set = x0_set;

figure;
feasibleSet = tube_mpc_diags.feasible_set;
infeasibleSet = tube_mpc_diags.infeasible_set;
scatter(x0_set(:,1), x0_set(:,2));
hold on
if ~isempty(feasibleSet)
    scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
end
title('tube MPC');
xlabel('x_1'); ylabel('x_2'); 

%% grid SLS MPC evaluation

% hyperparameter bounds
range = struct;
range.gamma.lb = 0.01; range.gamma.ub = 100;
range.beta.lb = 0.01; range.beta.ub = 100;
range.tau.lb = 0.01; range.tau.ub = 100;
bisectOptions = struct;
bisectOptions.tol = 1e-2;
[bounds, isFeasible] = mpc.BisectParams(range, bisectOptions);

feasibility_opt = struct;

gridDim = struct;
% grid search resolution 3 by 3 by 3
gridDim.num_beta = 5;
gridDim.num_gamma = 5;
gridDim.num_tau = 5;

gridparams_opt = struct;
gridparams_opt.verbose = 0;

% grid search range
range = struct;
range.gamma.lb = 0.01; range.gamma.ub = 100;
range.beta.lb = 0.01; range.beta.ub = 100;
range.tau.lb = 0.01; range.tau.ub = 100;

bisect_opt = struct;
bisect_opt.init.tau.lb = bounds.tau.lb;
bisect_opt.init.beta.lb = bounds.beta.lb;
bisect_opt.tol = 1e-2;
bisect_opt.verbose = 0;

feasibility_opt.range = range;
feasibility_opt.gridparams_opt = gridparams_opt;
feasibility_opt.gridDim = gridDim;
feasibility_opt.bisect_opt = bisect_opt;

start = tic;
[solCell, sets] = VerifyFeasibility(mpc, x0_set, feasibility_opt);
naive_sls_mpc_running_time = toc(start);

naive_sls_mpc_diags = struct;
naive_sls_mpc_diags.feasible_set = sets.feasibleSet;
naive_sls_mpc_diags.infeasible_set = sets.infeasibleSet;
naive_sls_mpc_diags.unverified_set = sets.unverifiedSet;
naive_sls_mpc_diags.running_time = naive_sls_mpc_running_time;
naive_sls_mpc_diags.x0_set = x0_set;

figure;
feasibleSet = naive_sls_mpc_diags.feasible_set;
infeasibleSet = naive_sls_mpc_diags.infeasible_set;
unverifiedSet = naive_sls_mpc_diags.unverified_set;
scatter(x0_set(:,1), x0_set(:,2));
hold on
if ~isempty(feasibleSet)
    scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
end
if ~isempty(unverifiedSet)
    scatter(unverifiedSet(:,1), unverifiedSet(:,2),'g*');
end
title('tube MPC');
xlabel('x_1'); ylabel('x_2'); 

%% save data
% load the following dataset to recover the Figure 1 in the paper
% save('data\conservatism_0723.mat');

%% plot 
Xc = mpc.state_constr;
Uc = mpc.input_constr;
XT = mpc.terminal_constr;
figure; hold on
Xc.plot('wire', true, 'linewidth', 2, 'linestyle', '-');
% Xc.plot('color', 'lightgreen', 'linestyle', 'none');
XT.plot('color', 'lightblue', 'linestyle', 'none');
% scatter(x0_set(:,1), x0_set(:,2), 4 );

aug_sls_mpc_feasible_set = aug_df_mpc_diags.feasible_set;
ROA_aug_sls = Polyhedron(aug_sls_mpc_feasible_set);
unif_df_mpc_feasible_set = unif_df_mpc_diags.feasible_set;
ROA_unif_df = Polyhedron(unif_df_mpc_feasible_set);
grid_sls_mpc_feasible_set = naive_sls_mpc_diags.feasible_set;
ROA_grid_sls = Polyhedron(grid_sls_mpc_feasible_set);
tube_mpc_feasible_set = tube_mpc_diags.feasible_set;
ROA_tube = Polyhedron(tube_mpc_feasible_set);

h1 = ROA_aug_sls.plot('wire', true, 'linewidth', 2, 'linestyle', '-', 'edgecolor', [0 0.4470 0.7410]);
h2 = ROA_unif_df.plot('wire', true, 'linewidth', 2, 'linestyle', '-.', 'edgecolor', 	[0.8500 0.3250 0.0980]);
h3 = ROA_grid_sls.plot('wire', true, 'linewidth', 2, 'linestyle', '--', 'edgecolor', [0.9290 0.6940 0.1250]);
h4 = ROA_tube.plot('wire', true, 'linewidth', 2, 'linestyle', ':', 'edgecolor', [0.4940 0.1840 0.5560]);
legend([h1, h2, h3, h4], {'aug-SLS', 'unif-df', 'grid-SLS', 'tube'}, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 18);

