function [result] = MPC_feasibility_compare(mpc, x0_set, Z)
%MPC_FEASIBILITY_COMPARE Summary of this function goes here
%   x0_set: x0 on each row

methods = {'aug_sls', 'unif_df', 'tube', 'grid_sls'};

result = struct;

aug_sls_diags = MPC_feasibility_check(mpc, x0_set, 'aug_sls');
result.aug_sls = aug_sls_diags;

unif_df_diags = MPC_feasibility_check(mpc, x0_set, 'unif_df');
result.unif_df = unif_df_diags;

tube_diags = MPC_feasibility_check(mpc, x0_set, 'tube', Z);
result.tube = tube_diags;

grid_sls_diags = MPC_feasibility_check(mpc, x0_set, 'grid_sls');
result.grid_sls = grid_sls_diags;

result.total_num = size(x0_set, 1);
result.methods = methods;
num_feasible_points = [];
for ii = 1:4
   method = methods{ii};
   num = size(result.(method).feasible_set,1);
   num_feasible_points = [num_feasible_points num]; 
end
result.num_feasible_points = num_feasible_points;
result.coverage = num_feasible_points/result.total_num;
end

