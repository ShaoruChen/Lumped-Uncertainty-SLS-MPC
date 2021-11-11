clear;
%% initialization

nx = 2;
nu = 1;

A = [1 0.15; 0.1 1];
B = [0.1; 1.1];

x0 = [-3; 4];

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

E_bar = [1 0; 0 1; -1 0; 0 -1];
e_bar = [8; 8; 8; 8];
Xc = Polyhedron(E_bar, e_bar);

% stage cost weights
Q = 10*eye(nx); R = eye(nu);

% uncomment to generate maximum robust positive invariant set
% opts = struct;
% opts.robust = 1; opts.minVol = 0.1;
% [RIS, diagnostic] = uncertain_system.robustInvariantSet(Xc, Uc, 20, opts);
% save('RIS', 'RIS');

% uncomment to generate the robust positive invariant set w.r.t. a local controller u = Kx
% opts = struct;
% opts.robust = 1; opts.minVol = 0.5; opts.plot = 0;

% [K, P] = uncertain_system.find_K_LQR(Q, R);
% [RIS, diagnostic] = uncertain_system.robustInvariantSetClosedLoop(Xc, Uc, 30, opts);

RIS_data = load('RIS.mat');
terminal_set = RIS_data.RIS;


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
uncertain_system.find_K_LQR(Q,R);
[Z_inv, isConverge] = uncertain_system.minInvSet(20);
save('mpc_data', 'mpc', 'Z_inv');


