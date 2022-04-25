function [diags] = MPC_feasibility_check(mpc, x0_set, method, Z)
%MPC_FEASIBILITY_CHECK Summary of this function goes here
% method: 'tube': tube MPC; 'aug_sls', 'grid_sls', 'unif_df'
% Z: the tube, provided when method = 'tube'
if nargin < 4
   if strcmp(method, 'tube')
      error('Tube has to be provided for the tube MPC method.');
   end
end

if strcmp(method, 'aug_sls')
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

    diags = struct;
    diags.feasible_set = feasibleSet;
    diags.infeasible_set = infeasibleSet;
    diags.running_time = aug_df_MPC_running_time;
    diags.x0_set = x0_set;
elseif strcmp(method,'unif_df')
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

    diags = struct;
    diags.feasible_set = feasibleSet;
    diags.infeasible_set = infeasibleSet;
    diags.running_time = uniform_df_MPC_running_time;
    diags.x0_set = x0_set;
elseif strcmp(method, 'tube')
    num_points = size(x0_set, 1);
    feasibleSet =[]; infeasibleSet = [];
    start_time = tic;

    verbose = 0;
    [tube_mpc_optimizer] = mpc.SolveTubeMPC(Z, 'optimizer', verbose);

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
        diags = struct;
        diags.feasible_set = feasibleSet;
        diags.infeasible_set = infeasibleSet;
        diags.running_time = tube_MPC_running_time;
        diags.x0_set = x0_set;
elseif strcmp(method, 'grid_sls')
    % grid search parameters can be adjusted
    mpc.x0 = zeros(mpc.nx, 1);
    % hyperparameter bounds
    range = struct;
    range.gamma.lb = 0.001; range.gamma.ub = 50;
    range.beta.lb = 0.001; range.beta.ub = 50;
    range.tau.lb = 0.001; range.tau.ub = 50;
    bisectOptions = struct;
    bisectOptions.tol = 1e-2;
    [bounds, isFeasible] = mpc.BisectParams(range, bisectOptions);

    feasibility_opt = struct;

    gridDim = struct;
    % grid search resolution 5 by 5 by 5
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

    diags = struct;
    diags.feasible_set = sets.feasibleSet;
    diags.infeasible_set = sets.infeasibleSet;
    diags.unverified_set = sets.unverifiedSet;
    diags.running_time = naive_sls_mpc_running_time;
    diags.x0_set = x0_set;
else
    error('MPC method not supported.');
end

end

