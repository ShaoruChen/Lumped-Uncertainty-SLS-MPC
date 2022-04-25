% Compare the feasible region of robust MPC methods on randomly generated systems.
% Run section postprocessing directly to recover the Figure 3 in the paper by 
% loading existing dataset.

clear;

%% initialization

nx = 2;
nu = 1;

total_num_expr = 20;
expr_sol = cell(1, total_num_expr);

for kk = progress(1:total_num_expr)
expr_result = struct;

A = (rand(nx, nx) - 0.5)*4;
B = (rand(nx, nu) - 0.5)*2;

expr_result.A = A; expr_result.B = B;


x0 = zeros(nx, 1);

eps_A = 0.1;
eps_B = 0.1;
sigma_w = 0.1;

system_params = struct;
system_params.A = A;
system_params.B = B;
system_params.x0 = x0;
system_params.eps_A = eps_A;
system_params.eps_B = eps_B;
system_params.sigma_w = sigma_w;

uncertain_system = UncertainLTISystem(system_params);

%% state and input constraints
Uc_vertices = [-4; 4];
Uc = Polyhedron(Uc_vertices);

E = [eye(nx); -eye(nx)]; e = [8*ones(nx, 1); 8*ones(nx, 1)];
Xc = Polyhedron(E, e);

Q = eye(nx); R = eye(nu);

terminal_set = [];

%% construct SLS MPC problem

MPC_data = struct;
MPC_data.uncertain_system = uncertain_system;

horizon = 5; 
MPC_data.horizon = horizon;
MPC_data.eps_A = eps_A; MPC_data.eps_B = eps_B; MPC_data.sigma_w = sigma_w;

MPC_data.Q = Q; MPC_data.R = R;
MPC_data.state_constr = Xc; 
MPC_data.input_constr = Uc;
MPC_data.terminal_constr = terminal_set;

mpc = SLSMPC(uncertain_system, MPC_data);

% find tubes for tube mpc
uncertain_system.find_K_LQR(Q, R);
[Z_inv, isConverge] = uncertain_system.minInvSet(50);
if isConverge ~= 1
   warning('Z_inv not converged'); 
end
expr_result.Z_inv = Z_inv; expr_result.isConverge = isConverge;

save('mpc_data_random', 'mpc', 'Z_inv');

%% 
X_grid = mpc.state_constr;
x0_set = X_grid.grid(10);

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

expr_result.aug_df_mpc_diags = aug_df_mpc_diags;
% 
% figure;
% feasibleSet = aug_df_mpc_diags.feasible_set;
% infeasibleSet = aug_df_mpc_diags.infeasible_set;
% scatter(x0_set(:,1), x0_set(:,2));
% hold on
% if ~isempty(feasibleSet)
%     scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
% end
% title('augmented disturbance feedback');
% xlabel('x_1'); ylabel('x_2'); 

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

expr_result.unif_df_mpc_diags = unif_df_mpc_diags;

% figure;
% feasibleSet = unif_df_mpc_diags.feasible_set;
% infeasibleSet = unif_df_mpc_diags.infeasible_set;
% scatter(x0_set(:,1), x0_set(:,2));
% hold on
% if ~isempty(feasibleSet)
%     scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
% end
% title('uniform disturbance feedback');
% xlabel('x_1'); ylabel('x_2'); 

%% tube MPC evaluation
num_points = size(x0_set, 1);
feasibleSet =[]; infeasibleSet = [];
start_time = tic;

verbose = 0;
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

expr_result.tube_mpc_diags = tube_mpc_diags;

% figure;
% feasibleSet = tube_mpc_diags.feasible_set;
% infeasibleSet = tube_mpc_diags.infeasible_set;
% scatter(x0_set(:,1), x0_set(:,2));
% hold on
% if ~isempty(feasibleSet)
%     scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
% end
% title('tube MPC');
% xlabel('x_1'); ylabel('x_2'); 

%% grid SLS MPC evaluation

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

expr_result.naive_sls_mpc_diags = naive_sls_mpc_diags;
expr_sol{kk} = expr_result;
save('random_example_search_temp_0730.mat');
% figure;
% feasibleSet = naive_sls_mpc_diags.feasible_set;
% infeasibleSet = naive_sls_mpc_diags.infeasible_set;
% unverifiedSet = naive_sls_mpc_diags.unverified_set;
% scatter(x0_set(:,1), x0_set(:,2));
% hold on
% if ~isempty(feasibleSet)
%     scatter(feasibleSet(:,1), feasibleSet(:,2),'r*');
% end
% if ~isempty(unverifiedSet)
%     scatter(unverifiedSet(:,1), unverifiedSet(:,2),'g*');
% end
% title('tube MPC');
% xlabel('x_1'); ylabel('x_2'); 

end

% save('random_example_search_0730.mat');

%% post processing
load('random_example_search_0730.mat');
aug_df_num = []; unif_df_num = []; tube_num = []; naive_sls_num = [];
for ii= 1:total_num_expr
    aug_df_num = [aug_df_num size(expr_sol{ii}.aug_df_mpc_diags.feasible_set,1)];
    unif_df_num = [unif_df_num size(expr_sol{ii}.unif_df_mpc_diags.feasible_set,1)];
    tube_num = [tube_num size(expr_sol{ii}.tube_mpc_diags.feasible_set,1)];
    naive_sls_num = [naive_sls_num size(expr_sol{ii}.naive_sls_mpc_diags.feasible_set,1)];
end

% figure; hold on
% plot(aug_df_num, '-');
% plot(unif_df_num, '-.');
% plot(tube_num, ':');
% plot(naive_sls_num, '--');
% 
% [~, ind] = sort(aug_df_num);
% figure; hold on
% plot(aug_df_num(ind), '-');
% plot(unif_df_num(ind), '-.');
% plot(tube_num(ind), ':');
% plot(naive_sls_num(ind), '--');
% xlabel('example number'); ylabel('coverage');

%% combine two data sets
data = load('random_example_search_0729.mat');
expr_sol_1 = data.expr_sol;
total_num = data.total_num_expr;

for ii= 1:total_num
    aug_df_num = [aug_df_num size(expr_sol_1{ii}.aug_df_mpc_diags.feasible_set,1)];
    unif_df_num = [unif_df_num size(expr_sol_1{ii}.unif_df_mpc_diags.feasible_set,1)];
    tube_num = [tube_num size(expr_sol_1{ii}.tube_mpc_diags.feasible_set,1)];
    naive_sls_num = [naive_sls_num size(expr_sol_1{ii}.naive_sls_mpc_diags.feasible_set,1)];
end


[~, ind] = sort(aug_df_num);
figure; hold on
plot(aug_df_num(ind)/100, '-', 'LineWidth',2 );
plot(unif_df_num(ind)/100, '-.', 'LineWidth',2 );
plot(tube_num(ind)/100, ':', 'LineWidth',2 );
plot(naive_sls_num(ind)/100, '--', 'LineWidth',2 );
legend('aug-SLS', 'unif-df', 'tube', 'grid-SLS', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('example number', 'Interpreter', 'latex', 'FontSize', 18); 
ylabel('coverage', 'Interpreter', 'latex', 'FontSize', 18);

is_stable = [];
for ii = 1:total_num_expr
   is_stable = [is_stable max(eig(expr_sol{ii}.A))< 1];
end

for ii = 1:total_num
       is_stable = [is_stable max(eig(expr_sol_1{ii}.A))< 1];
end
    
