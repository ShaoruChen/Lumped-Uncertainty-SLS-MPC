function [mpc, Z_tube] = problem_generator(eps_A, eps_B, sigma_w, horizon)
%PROBLEM_GENERATOR generate robust MPC problem instances

nx = 2;
nu = 1;

A = [1 0.15; 0.1 1];
B = [0.1; 1.1];

% construct polytopic model uncertainty
Delta_A_1 = [eps_A 0; 0 0]; Delta_A_2 = [-eps_A 0; 0 0]; 
Delta_A_vertices = {Delta_A_1, Delta_A_2};

Delta_B_1 = [0; -eps_B]; Delta_B_2 = [0; eps_B];
Delta_B_vertices = {Delta_B_1, Delta_B_2};

x0 = [-3; 4];

system_params = struct;
system_params.A = A;
system_params.B = B;
system_params.x0 = x0;
system_params.eps_A = eps_A;
system_params.eps_B = eps_B;
system_params.sigma_w = sigma_w;

system_params.Delta_A_vertices = Delta_A_vertices;
system_params.Delta_B_vertices = Delta_B_vertices;
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
opts = struct;
opts.robust = 1; opts.minVol = 0.1;
[RIS, diagnostic] = uncertain_system.robustInvariantSet(Xc, Uc, 100, opts);
if diagnostic.converge == 0
   warning('Maximum robust invariant set iterations fail to converge. '); 
end

terminal_set = RIS;

% uncomment to generate the robust positive invariant set w.r.t. a local controller u = Kx
% opts = struct;
% opts.robust = 1; opts.minVol = 0.5; opts.plot = 0;

% [K, P] = uncertain_system.find_K_LQR(Q, R);
% [RIS, diagnostic] = uncertain_system.robustInvariantSetClosedLoop(Xc, Uc, 30, opts);


%% construct SLS MPC problem

MPC_data = struct;
MPC_data.uncertain_system = uncertain_system;

MPC_data.horizon = horizon;
MPC_data.eps_A = eps_A; MPC_data.eps_B = eps_B; MPC_data.sigma_w = sigma_w;

MPC_data.Q = Q; MPC_data.R = R;
MPC_data.state_constr = Xc; 
MPC_data.input_constr = Uc;
MPC_data.terminal_constr = terminal_set;

MPC_data.Delta_A_vertices = Delta_A_vertices;
MPC_data.Delta_B_vertices = Delta_B_vertices;

mpc = SLSMPC(uncertain_system, MPC_data);

% find tubes for tube mpc
uncertain_system.find_K_LQR(Q,R);
[Z_tube, isConverge] = uncertain_system.minInvSet(100);
if isConverge == 0
    warning('Disturbance invariant set iterations fail to converge.');
end

end

