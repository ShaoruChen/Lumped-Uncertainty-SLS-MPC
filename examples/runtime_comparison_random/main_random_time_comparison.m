clear;
%% initialization

nx = 2;
nu = 1;

total_num_expr = 100;
expr_sol = cell(1, total_num_expr);
for kk = progress(1:total_num_expr)
    
expr_result = struct;

M = 4;
A = (rand(nx, nx) - 0.5)*M;
B = (rand(nx, nu) - 0.5)*2;

expr_result.A = A; expr_result.B = B;

% sys_data = load('sys_dyn.mat');
% A = sys_data.A; B = sys_data.B;

x0 = zeros(nx, 1);

eps_A = 0.05;
eps_B = 0.05;
sigma_w = 0.05;

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

horizon = 10; 
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

x0 = [2 -1]';
mpc.x0 = x0;

%% augmented disturbance feedback MPC evaluation
opt = struct;
opt.solver = 'mosek'; opt.verbose = 0;
[aug_df_sol] = mpc.SolveAugDistFeedbackSLSMPC('value', opt);

aug_df_mpc_diags = struct;
aug_df_mpc_diags.sol = aug_df_sol
aug_df_mpc_diags.solver_time = aug_df_sol.solver_time;
expr_result.aug_df_mpc_diags = aug_df_mpc_diags;


%% uniform disturbance feedback MPC evaluation
verbose = 0;
[unif_df_sol] = mpc.SolveUniformDistFeedbackMPC('value', verbose);
unif_df_mpc_diags = struct;
unif_df_mpc_diags.sol = unif_df_sol;
if ~isfield(unif_df_sol, 'solver_time')
    unif_df_sol.solver_time = nan;
end
unif_df_mpc_diags.solver_time = unif_df_sol.solver_time;
expr_result.unif_df_mpc_diags = unif_df_mpc_diags;

%% tube MPC evaluation

verbose = 0;
[tube_mpc_sol] = mpc.SolveTubeMPC(Z_inv, 'value', verbose);

tube_mpc_diags = struct;
tube_mpc_diags.sol = tube_mpc_sol;
tube_mpc_diags.solver_time = tube_mpc_sol.solver_time;
expr_result.tube_mpc_diags = tube_mpc_diags;

%% grid SLS MPC evaluation
[naive_sls_sol] = mpc.SolveSLSMPCAuto()

naive_sls_diags = struct;
naive_sls_diags.sol = naive_sls_sol;
naive_sls_diags.solver_time = naive_sls_sol.solver_time;
expr_result.naive_sls_diags = naive_sls_diags;

expr_sol{kk} = expr_result;

save('temp_data.mat');

end

save('random_exampl_time_comparison_0906.mat');

%% post processing
aug_df_count = zeros(1, total_num_expr);
unif_df_count = zeros(1, total_num_expr);
tube_count = zeros(1, total_num_expr);
naive_sls_count = zeros(1, total_num_expr);

for ii = 1:total_num_expr
    expr_result = expr_sol{ii};
    
    % aug_df_mpc 
    aug_df_diags = expr_result.aug_df_mpc_diags;
    if aug_df_diags.sol.status == 0 
        aug_df_count(ii) = aug_df_diags.solver_time;
    else
        aug_df_count(ii) = nan;
    end
    
    % unif_df_mpc
    unif_df_diags = expr_result.unif_df_mpc_diags;
    if unif_df_diags.sol.status == 0
        unif_df_count(ii) = unif_df_diags.solver_time;
    else
        unif_df_count(ii) = nan;
    end
    
    tube_diags = expr_result.tube_mpc_diags;
    if tube_diags.sol.status == 0 
        tube_count(ii) = tube_diags.solver_time;
    else
        tube_count(ii) = nan;
    end
    
    naive_sls_diags = expr_result.naive_sls_diags;
    if naive_sls_diags.sol.status == 0
        naive_sls_count(ii) = naive_sls_diags.solver_time;
    else
       naive_sls_count(ii) = nan; 
    end
    
end

effective_aug_df_count = aug_df_count(~isnan(aug_df_count));
effective_unif_df_count = unif_df_count(~isnan(unif_df_count));
effective_tube_count = tube_count(~isnan(tube_count));
effective_naive_sls_count = naive_sls_count(~isnan(naive_sls_count));

data = [effective_aug_df_count effective_unif_df_count effective_tube_count effective_naive_sls_count];
group = [zeros(size(effective_aug_df_count)) ones(size(effective_unif_df_count)) ...
       2*ones(size(effective_tube_count))  3*ones(size(effective_naive_sls_count))];

figure;
boxplot(data', group',  'labels', {'aug-SLS-MPC', 'unif-df-MPC', 'tube-MPC', 'grid-SLS-MPC'});
ax = gca;
ax.YAxis.Scale ="log";

set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'aug-SLS-MPC', 'unif-df-MPC', 'tube-MPC', 'grid-SLS-MPC'},'FontSize',12);
ylabel('solver time [sec]', 'Interpreter', 'Latex', 'FontSize', 18);


%% count unstale systems
eig_max = zeros(1, total_num_expr);
for ii = 1:total_num_expr
    expr_result = expr_sol{ii};
    A_mat = expr_result.A;
    eig_max(ii) = max(abs(eig(A_mat)));
end
sum(eig_max > 1)
