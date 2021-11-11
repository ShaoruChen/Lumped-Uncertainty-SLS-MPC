clear;

%% test different robust MPC methods

data = load('mpc_data.mat');
mpc = data.mpc;
Z_inv = data.Z_inv;

%% augmented disturbance feedback MPC
[aug_df_sol] = mpc.SolveAugDistFeedbackSLSMPC();

[traj_set] = mpc.SimulateTrajwithSystemResponses(20);
PlotFcns.plot_SLS_MPC_traj(mpc, traj_set, 'SLS MPC');

%% uniform disturbance feedback MPC
[unif_df_sol] = mpc.SolveUniformDistFeedbackMPC(2);

num_traj = 20; proj_dim = [1 2];
[unif_df_traj_set] = mpc.SimulateTrajwithUniformDistFeedbackMPC(unif_df_sol, num_traj);
PlotFcns.plot_uniform_dist_feedback_MPC_traj(mpc, unif_df_traj_set, 'Unif. DF MPC', proj_dim);
%% tube MPC test
% find tubde

[tube_sol] = mpc.SolveTubeMPC(Z_inv, 'value', 2);

num_traj = 5;
proj_dim = [1 2];
[tube_traj_set] = mpc.SimulateTrajwithTubeMPC(tube_sol, num_traj)
PlotFcns.plot_tube_MPC_traj(mpc, tube_sol, tube_traj_set,'Tube MPC', proj_dim);

%% SLS MPC
[naive_sls_sol] = mpc.SolveSLSMPCAuto();
num_traj = 5; proj_dim = [1 2];
[naive_sls_traj_set] = mpc.SimulateTrajwithSystemResponses(num_traj);

PlotFcns.plot_SLS_MPC_traj(mpc, naive_sls_traj_set, 'SLS MPC', proj_dim);
