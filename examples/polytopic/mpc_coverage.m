function [mpc_diags] = mpc_coverage(mpc, x0_set, type, params)
%MPC_COVERAGE evaluate the feasible domain coverage of robust MPC methods

if nargin < 4
    params = [];
end

num_points = size(x0_set, 1);
feasibleSet =[]; infeasibleSet = [];

if strcmp(type, 'polytopic_SLS_MPC')
opt = struct;
opt.solver = 'mosek'; opt.verbose = 0; 
opt.norm_type = inf;
[poly_sls_mpc_optimizer] = mpc.SolvePolytopicSLSMPC('optimizer', opt);

start_time = tic;
for ii = progress(1:num_points)
    init_x = x0_set(ii, :)';
    [sol_value, errorcode] = poly_sls_mpc_optimizer(init_x);
    
    if errorcode == 0
        feasibleSet = [feasibleSet; init_x'];
    else
        infeasibleSet = [infeasibleSet; init_x'];
    end
end
poly_sls_mpc_running_time = toc(start_time);

mpc_diags = struct;
mpc_diags.feasible_set = feasibleSet;
mpc_diags.infeasible_set = infeasibleSet;
mpc_diags.running_time = poly_sls_mpc_running_time;
mpc_diags.x0_set = x0_set;
mpc_diags.coverage = size(feasibleSet, 1)/num_points;
mpc_diags.mpc = mpc;

elseif strcmp(type, 'tube_MPC')
    % For tube MPC, params should give a template tube cross section.
    Z_tube = params;
    [tube_mpc_optimizer] = mpc.SolvePolytopicTubeMPC(Z_tube, 'optimizer');
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

    mpc_diags = struct;
    mpc_diags.feasible_set = feasibleSet;
    mpc_diags.infeasible_set = infeasibleSet;
    mpc_diags.running_time = tube_MPC_running_time;
    mpc_diags.x0_set = x0_set;
    mpc_diags.Z = Z_tube;
    mpc_diags.coverage = size(feasibleSet,1)/num_points;
    mpc_diags.mpc = mpc;
else
    error('The specified method not supported.');
end

end

