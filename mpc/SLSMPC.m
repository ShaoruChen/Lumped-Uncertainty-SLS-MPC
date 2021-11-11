classdef SLSMPC < handle
    %SLSMPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        horizon;
        uncertain_system;
        nx; nu; na; nb;
        eps_A; eps_B; sigma_w;
        x0; 
        x_eq; % goal state
        
        Q; R; Q_T;
        state_constr; input_constr; 
        terminal_constr = [];           
        Phi_x = []; Phi_u = [];
        
        sigma_seq = [];
    end
    
    methods (Access = public)
        %% initialize the SLSMPC problem
        function obj = SLSMPC(uncertain_system, params)
            obj.uncertain_system = uncertain_system;
            obj.nx = uncertain_system.nx;
            obj.nu = uncertain_system.nu;
            obj.x0 = uncertain_system.x0;
            
            obj.horizon = params.horizon;
            obj.eps_A = params.eps_A;
            obj.eps_B = params.eps_B;
            obj.sigma_w = params.sigma_w;
            
            obj.Q = params.Q;
            obj.R = params.R;
            % set the equilibrium point
            if isfield(params, 'Q_T')
                obj.Q_T = params.Q_T;
            else
                obj.Q_T = params.Q;
            end
            
            if isfield(params, 'x_ref')
                obj.x_eq = params.x_ref;
            else
                obj.x_eq = zeros(obj.nx, 1);
            end
            
            obj.state_constr = params.state_constr;
            obj.input_constr = params.input_constr;
            if isfield(params,'terminal_constr') && ~isempty(params.terminal_constr)
                obj.terminal_constr = params.terminal_constr;
            else
                obj.terminal_constr = params.state_constr;
            end

            [unif_sigma] = obj.Find_uniform_ub_aug_dist();
            obj.sigma_seq = unif_sigma*ones(1, obj.horizon);
        end
        
        function set_x0(obj, x0)
            % reset the initial state of the MPC problem
            obj.x0 = x0;
            obj.uncertain_system.x0 = x0;
        end
        
        function [unif_sigma] = Find_uniform_ub_aug_dist(obj)
            % find a uniform upper bound on the augmented disturbances over
            % the state space and control input space
            sigma_w = obj.sigma_w; eps_A = obj.eps_A; eps_B = obj.eps_B;
            
            max_x_inf_norm = max(abs(obj.state_constr.V(:)));
            max_u_inf_norm = max(abs(obj.input_constr.V(:)));
            unif_sigma = eps_A*max_x_inf_norm + eps_B*max_u_inf_norm + sigma_w;
            
        end
        
        %% solve the SLS convex tightening of predictive control
        function [sol] = SolveSLSMPC(obj, hyperparams, opt) 
            % solve SLS MPC with given hyperparameters
            % Inputs:
            %   params: struct of hyperparameters. Fields: .beta, .gamma,
            %   .tau, .alpha (optional)
            %   opt: struct, fields: .verbose, .keepGamma, .keepAffine,
            %   .keepTau, .keepRobustPerf, .keepSafety, .keepObj
            
            if nargin < 3
                % options to choose which constraint to mute
                opt.verbose = 2;
                opt.keepGamma = 1;
                opt.keepAffine = 1;
                opt.keepTau = 1;
                opt.keepRobustPerf = 1;
                opt.keepSafety = 1;
                opt.keepObj = 1;
            end
             
            tau = hyperparams.tau; gamma = hyperparams.gamma; beta = hyperparams.beta;
            
            n = obj.nx; p = obj.nu; na = obj.na; nb = obj.nb;
            T = obj.horizon;
            
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            eps_sum = eps_A + eps_B;
            
            uncertain_system = obj.uncertain_system;
            x0 = obj.x0;

            sigma_w_aug = sigma_w;
            
            % construct system response variables 
			Phi_x = sdpvar( (T + 1) * n, (T + 1) * n, 'full');
			Phi_u = sdpvar( (T + 1) * p, (T + 1) * n, 'full');
            % extract matrix blocks
            Phi = [Phi_x; Phi_u];
            Phi_x_0 = Phi_x(:,1:n);
            Phi_u_0 = Phi_u(:,1:n);
            Phi_x_w = Phi_x(:,n+1:end);
            Phi_u_w = Phi_u(:,n+1:end);
            Phi_0 = Phi(:,1:n);
            Phi_w = Phi(:,n+1:end);
            
            % Construct the objective function using nominal cost
			sqrtQ = sqrtm(obj.Q);
			sqrtR = sqrtm(obj.R);
            sqrtQ_T = sqrtm(obj.Q_T);

			Qset = cell(T + 1, 1); Rset = cell(T + 1, 1);
            % no penalty on the initial state
			Qset{1} = zeros(size(sqrtQ));
			for i = 2 : T 
				Qset{i} = sqrtQ;
			end
			Qset{T + 1} = sqrtQ_T;

			for i = 1:T
				Rset{i} = sqrtR;
            end
            
            % attention: no penalty on u_T
			Rset{T + 1} = zeros(size(sqrtR));
%             Rset{T + 1} = eye(size(sqrtR));

			cost = 0;
            if opt.keepObj
                for i  = 1 : T+1
                    % attention: add reference target state
                    cost = cost + norm(Qset{i}*(Phi_x((i - 1)*n + 1: i*n, 1:n )*x0 ), 2)^2 ...
                           + norm(Rset{i}*(Phi_u((i - 1)*p + 1: i*p, 1:n )*x0), 2)^2;      
                end
            end
            
            % add constraints to the problem
            constr = [];
            
            % Structure constraint of Phi_x and Phi_u
			for k = 1 : T
				constr = [constr, Phi_x( (k - 1)*n + 1: k*n, k*n + 1: end) == zeros(n, (T + 1 - k)*n)];
            end

			for k = 1 : T
				constr = [constr, Phi_u( (k - 1)*p + 1: k*p, k*n + 1: end) == zeros(p, (T + 1 - k)*n)];
            end
            
            % define the block downshift operator
            Z = kron(diag(ones(1,T),-1), eye(n));
            A_block = blkdiag(kron(eye(T), uncertain_system.A), zeros(n, n));
            B_block = blkdiag(kron(eye(T), uncertain_system.B), zeros(n, p));
            ZA_block = Z*A_block; ZB_block = Z*B_block;
			
			Id = eye((T + 1)*n);
            
            if opt.keepAffine
                constr = [constr, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Id];
            end
            
            % state and input constraints
            % % construct F
			if ~isempty(obj.state_constr)
				Fx = obj.state_constr.A; bx = obj.state_constr.b;
				nFx = size(Fx, 1); nbx = length(bx);
            else
                warning('state constraints are required.');
            end
            
            if ~isempty(obj.terminal_constr)
				Ft = obj.terminal_constr.A; bt = obj.terminal_constr.b;
				nFt = size(Ft, 1); nbt = length(bt);
            else 
                Ft = Fx; bt = bx; nFt = nFx; nbt = nbx;
            end
            
            if ~isempty(obj.input_constr)
				Fu = obj.input_constr.A; bu = obj.input_constr.b;
				nFu = size(Fu, 1); nbu = length(bu);
            else
                warning('must have input constraints');
            end
            
            % state constraint
            F = blkdiag(Fx);
            for j = 1:T-1
                F = blkdiag(F, Fx);
            end
            % terminal state constraint
            F = blkdiag(F, Ft);
            % control input constraint
            for j = 1:T
                F = blkdiag(F, Fu);
            end
            % concatenate all input, state and terminal constraints.
            % constraint on u_T will not be included in the optimal control
            % problem
            F = blkdiag(F, zeros(nFu, p));
            b = [kron(ones(T,1), bx); kron(ones(1,1), bt); kron(ones(T,1),bu); kron(ones(1,1), zeros(size(bu)))];
            
            nFconstr = size(F,1);
            assert(nFconstr == T*nFx + nFt + (T+1)*nFu);
            assert(size(b,1) == T*nFx + nFt + (T+1)*nFu);
            assert(size(b,2) == 1);    
            
            
            for j = 1:nFconstr - nFu % no constraint on u_T
                if opt.keepSafety
                   % since F(j,:)*Phi_w is a row vector, 1-norm is used
                   % here which corresponds to the ell_inf induced norm for
                   % matrices
                   constr = [constr, F(j,:)*Phi_0*x0 + norm(F(j,:)*Phi_w, 1)*(1 - tau^T)/(1 - tau)*gamma + beta*sigma_w_aug  <= b(j)];
                end

                if sigma_w_aug ~= 0
                    if opt.keepRobustPerf
                        constr = [constr, norm(F(j,:)*Phi_w, 1) + beta*norm(eps_sum*Phi_w, inf) <= beta];
                    end
                end
            end
            
            if opt.keepTau
                Phi_temp_1 = 2*[eps_A*Phi_x_w; eps_B*Phi_u_w];
                constr = [constr, norm(Phi_temp_1, inf) <= tau];
            end
            
            if opt.keepGamma
                Phi_temp_2 = 2*[eps_A*Phi_x_0; eps_B*Phi_u_0];
                constr = [constr, norm(Phi_temp_2*x0, inf) <= gamma];
            end
            
             % solve the problem
             if isfield(opt, 'solver')
                 solver = opt.solver;
             else
                 solver = 'mosek';
             end
             
             if isfield(opt, 'verbose')
                 verbose = opt.verbose;
             else
                 verbose = 2;
             end
             
			ops = sdpsettings('verbose', verbose, 'solver', solver);
            solution = optimize(constr, cost, ops);      
      
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
			end

			sol = struct;
			Phi_x_val = value(Phi_x); Phi_u_val = value(Phi_u);
            fval = value(cost);
                
            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            sol.fval = fval; 
            sol.tau = tau; sol.gamma = gamma; sol.beta = beta;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;

            obj.Phi_x = Phi_x_val; obj.Phi_u = Phi_u_val;  
            yalmip('clear');% clear all yalmip variables

        end

         %% bisect to find lower and upper bounds of the hyperparameters gamma, beta, tau
        function [bounds, isFeasible, solvertime] = BisectParams(obj, range, options)
            % use bisection to find lower and upper bounds of the
            % hyperparameters used in SLS MPC
            % Inputs:
            %   range: struct, for the form range.*.lb and range.*.ub where
            %          * = gamma, beta, tau
            %   options: struct, fields: options.init, options.tol,
            %   options.verbose. options.init contains given upper or lower
            %   bounds on the hyerparameters.
            if ~isfield(options, 'init')
                init = struct;
            else
                init = options.init;
            end
            
            if ~isfield(options, 'tol')
                tol = 5*1e-3;
            else 
                tol = options.tol;
            end
            
            if ~isfield(options, 'verbose')
                verbose = 0;
            else
                verbose = options.verbose;
            end
            
            % hyperparams for obj.SolveSLSMPC
            params = struct;
            params.beta = 1;
            params.gamma = 0.2;
            params.tau = 0.2;

            isFeasible = 1;            
            bounds = struct;
            solvertime = 0;
            
            % find lower bounds for beta, gamma, tau first
            % Find lower bound on gamma
            fprintf('Finding lower bound on gamma...\n');

            if isfield(init,'gamma') && isfield(init.gamma, 'lb') && ~isempty(init.gamma.lb)
                bounds.gamma.lb = init.gamma.lb;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 1;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.gamma.lb; ub = range.gamma.ub;
                mid = (lb + ub)/2;

                while ub - lb > tol
                    params.gamma = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.gamma.lb = mid;
            end

            % Find lower bound on beta
            fprintf('Finding lower bound on beta...\n');

            if isfield(init,'beta') && isfield(init.beta, 'lb') && ~isempty(init.beta.lb)
                bounds.beta.lb = init.beta.lb;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 1;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.beta.lb; ub = range.beta.ub;
                mid = (lb + ub)/2;
                while ub - lb > tol
                    params.beta = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.beta.lb = mid;
            end
            
            % Find lower bound on tau
            fprintf('Finding lower bound on tau...\n');

            if isfield(init,'tau') && isfield(init.tau, 'lb') && ~isempty(init.tau.lb)
                bounds.tau.lb = init.tau.lb;  
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 1;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 0;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;

                while ub - lb > tol
                    params.tau = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        lb = mid; 
                        mid = (lb + ub)/2;
                    else
                        ub = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.tau.lb = mid;      
            end
            
            % find upperbounds for beta, gamma, tau
            % Find upper bound on beta
            fprintf('Finding upper bound on beta...\n');

            if isfield(init,'beta') && isfield(init.beta, 'ub') && ~isempty(init.beta.ub)
                bounds.beta.ub = init.beta.ub;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 1;
                opt.keepObj = 0;
            
                lb = range.beta.lb; ub = range.beta.ub;
                mid = (lb + ub)/2;
                params.tau = bounds.tau.lb; params.gamma = bounds.gamma.lb;
                while ub - lb > tol
                    params.beta = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        ub = mid; 
                        mid = (lb + ub)/2;
                    else
                        lb = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.beta.ub = mid; 
            end
            
            if bounds.beta.ub < bounds.beta.lb
                isFeasible = 0;
                return
            end
                       
            % find upperbound on tau
            fprintf('Finding upper bound on tau...\n');

            if isfield(init,'tau') && isfield(init.tau, 'ub') && ~isempty(init.tau.ub)
                bounds.tau.ub = init.tau.ub;
            else
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 1;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;
                params.gamma = bounds.gamma.lb; params.beta = bounds.beta.lb;
                while ub - lb > tol
                    params.tau = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        ub = mid; 
                        mid = (lb + ub)/2;
                    else
                        lb = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.tau.ub = mid;   
            end
            
            if bounds.tau.ub < bounds.tau.lb
                isFeasible = 0;
                return
            end
            
            % find upperbound on gamma
            fprintf('Finding upper bound on gamma...\n');

            if isfield(init,'gamma') && isfield(init.gamma, 'ub') && ~isempty(init.gamma.ub)
                bounds.gamma.ub = init.gamma.ub;
            else
%                 tempup = (1 - bounds.tau.ub^obj.horizon)/(1 - bounds.tau.ub);
%                 templow = (1 - bounds.tau.lb^obj.horizon)/(1 - bounds.tau.lb);
%                 bounds.gamma.ub = tempup*bounds.gamma.lb/templow;
                opt = struct;
                opt.verbose = verbose;
                opt.keepGamma = 0;
                opt.keepAffine = 1;
                opt.keepTau = 0;
                opt.keepRobustPerf = 0;
                opt.keepSafety = 1;
                opt.keepObj = 0;

                lb = range.tau.lb; ub = range.tau.ub;
                mid = (lb + ub)/2;
                params.tau = bounds.tau.lb; params.beta = bounds.beta.lb;
                while ub - lb > tol
                    params.gamma = mid;
                    sol = obj.SolveSLSMPC(params, opt);
                    solvertime = solvertime + sol.solution.solvertime;
                    if sol.status ~= 0
                        ub = mid; 
                        mid = (lb + ub)/2;
                    else
                        lb = mid;
                        mid = (lb + ub)/2;
                    end    
                end
                bounds.gamma.ub = mid;   
                
            end
            
            if bounds.gamma.ub < bounds.gamma.lb
                isFeasible = 0;
                return
            end
            
            if (bounds.gamma.ub - bounds.gamma.lb < 1e-4) || (bounds.beta.ub - bounds.beta.lb < 1e-4) ...
                      || (bounds.tau.ub - bounds.tau.lb < 1e-4)
               isFeasible = 0;
               return
            end
            fprintf('Bisection finished.\n');
        end
        
          %% grid search for a feasible set of hyperparameters
        function [feasibleParams, solvertime] = GridSearchParams(obj, bounds, gridDim, options)
            % Grid search for a feasible set of hyperparameters (gamma,
            % beta, tau).
            % Inputs:
            %   bounds: struct, fields: bounds.*.lb, bounds.*.ub, with * = beta, tau,
            %   gamma.
            %   gridDim: struct, fields: gridDim.num_* with * = beta, tau, gamma.
            %   options: struct, fields: options.verbose
            
            if ~isfield(options, 'verbose')
                verbose = 0;
            else
                verbose = options.verbose;
            end
                   
            feasibleParams = [];
            solvertime = 0;
            
            solveOpt = struct;
            solveOpt.verbose = verbose;
            solveOpt.keepGamma = 1;
            solveOpt.keepAffine = 1;
            solveOpt.keepTau = 1;
            solveOpt.keepRobustPerf = 1;
            solveOpt.keepSafety = 1;
            solveOpt.keepObj = 0;
            
            num_tau = gridDim.num_tau; 
            num_gamma = gridDim.num_gamma; 
            num_beta = gridDim.num_beta;
            
            betaRange = linspace(bounds.beta.lb, bounds.beta.ub, num_beta + 2);
            gammaRange = linspace(bounds.gamma.lb, bounds.gamma.ub, num_gamma + 2);
            tauRange = linspace(bounds.tau.lb, bounds.tau.ub, num_tau + 2);
            
            % discard samples right on the boundary (they tend to be
            % infeasible)
            for beta = betaRange(2:end-1)
                for gamma = gammaRange(2:end-1)
                    for tau = tauRange(2:end-1)
                        params = struct;
                        params.beta = beta; params.gamma = gamma; params.tau = tau;
         
                        fprintf('grid search for SLS MPC feasible params: %f / %d, %f / %d, %f / %d \n', ...
                        find(betaRange == beta), length(betaRange)-2, find(gammaRange == gamma), length(gammaRange)-2, ...
                        find(tauRange == tau), length(tauRange)-2);
                        [sol] = obj.SolveSLSMPC(params, solveOpt); 
                        fprintf('solution code: %d \n', sol.status);
                        solvertime = solvertime + sol.solution.solvertime;
                        % return on finding the first feasible set
                        if sol.status == 0
                            feasibleParams = params;
                            fprintf('Feasible parameters found!\n');
                            return
                        end
                    end
                end
            end 
            fprintf('Grid search terminated. Good luck next time! \n');
        end
        
        %% solve SLS MPC pipepline
        function [sol] = SolveSLSMPCAuto(obj)
            % implement the whole pipeline of solving SLS MPC
            sol = struct;
            solver_time = 0;
            % search hyperparameter bounds
            range = struct;
            range.gamma.lb = 0.01; range.gamma.ub = 100;
            % fixme
            range.beta.lb = 0.01; range.beta.ub = 100;
            range.tau.lb = 0.01; range.tau.ub = 100;
            bisectOptions = struct;
            bisectOptions.tol = 1e-2;
            [bounds, isFeasible, solver_time_bisection] = obj.BisectParams(range, bisectOptions);
            sol.bounds = bounds;
            solver_time = solver_time + solver_time_bisection;
            
            if ~isFeasible
              sol.status = 1;
              sol.solver_time = solver_time;
              return
            end
            
            % grid search for a feasible hyperparameter
            gridDim = struct;
            gridDim.num_beta = 5;
            gridDim.num_gamma = 5;
            gridDim.num_tau = 5;
            options = struct;
            options.verbose = 0;
            [feasibleParams, solver_time_grid] = obj.GridSearchParams(bounds, gridDim, options)
            solver_time = solver_time + solver_time_grid;
            if isempty(feasibleParams)
                sol.status = 1;
                sol.solver_time = solver_time;
                return 
            end
            
            % solve SLS MPC
            [sol] = obj.SolveSLSMPC(feasibleParams);
            sol.bounds = bounds;
            sol.solver_time = solver_time + sol.solution.solvertime;
        end
        
         %% tube MPC approach
        function [sol] = SolveTubeMPC(obj, Z, type, verbose)
            % solve robust OCP through tube MPC
            % details see "Robust model predictive control using tubes" by
            % W.Langson, I.Chryssochoos, S.V.Rakovic, D.Q.Mayne
            % Z: a Polyhderon instance, the shape of the tube
            % type = 'value' to find the solution of tube MPC
            % type = 'optimizer' to generate a Yalmip optimizer with x_0 as
            % the input parameter 
            
            if nargin < 3
                type = 'value';
                verbose = 2;
            elseif nargin < 4
                verbose = 2;
            end
            
            fprintf('Solving tube MPC started... \n');

            uncertain_system = obj.uncertain_system;
            Ahat = uncertain_system.A; Bhat = uncertain_system.B;
            nx = obj.nx; nu = obj.nu;
            Q = obj.Q; R = obj.R; P = obj.Q_T;
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            T = obj.horizon;
            
            x0 = obj.x0; 
            Xc = obj.state_constr; 
            Uc = obj.input_constr;
            
            Wc = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
            Xf = obj.terminal_constr;
            
            Aw = Wc.A; bw = Wc.b;

            % shape of the tubes
            Z_vertices = Z.V;
            assert(size(Z_vertices,2) == nx);
            
            nJ = size(Z_vertices, 1);
            % set of vertices of Z
            Z_Vset = mat2cell(Z_vertices', [size(Z_vertices, 2)], ones(1, nJ));
            
            % find the H representation of Z
            Az = Z.A; bz = Z.b;
            assert(isempty(Z.Ae));
            
            % for computing the Minkowski difference
            nz = size(Az, 1);
            max_W = zeros(nz, 1);
            if sigma_w ~= 0
                opts = optimoptions(@linprog,'Display', 'off');
                for ii = 1:nz
                   [~, fval] = linprog(-Az(ii,:)', Aw, bw,[],[],[],[],opts );
                   max_W(ii) = -fval;
                end
            end
            
            % create optimization variables
            
           if strcmp(type, 'value')
                x0 = obj.x0;    
            elseif strcmp(type, 'optimizer')
                x0 = sdpvar(nx, 1);
           end
            
           if T > 1
                alpha = sdpvar(ones(1,T), ones(1,T));
           else
               alpha = {sdpvar(1,1)};
           end
                
            if T > 1
                zc = sdpvar(repmat(nx, 1, T), repmat(1, 1, T));
            else
                zc = {sdpvar(repmat(nx, 1, T), repmat(1, 1, T))};
            end
            u0 = sdpvar(nu, 1);
            for ii = 1:T
                % models {u_1, u_2, ..., u_T}
                U{ii} = sdpvar(repmat(nu, 1, nJ), repmat(1, 1, nJ));
                % X{i}{j}: the j-th vertex of the tube X_i
                % X_i contains x_i in the prediction
                X{ii} = sdpvar(repmat(nx, 1, nJ), repmat(1, 1, nJ));
            end

            for ii = 1:T
                for jj = 1:nJ
                   X{ii}{jj} = zc{ii} + alpha{ii}*Z_Vset{jj};
                end
            end
            
            % vertices of the inf-norm ball of Delta_A and Delta_B
            [Delta_A_vertices] = FindMatVertices(size(Ahat), eps_A);
            [Delta_B_vertices] = FindMatVertices(size(Bhat), eps_B);
            num_Delta_A = length(Delta_A_vertices);
            num_Delta_B = length(Delta_B_vertices);
            
            % construct the optimization problem
            fprintf('Constructing constraints ... \n');
            F = [ismember(u0, Uc)];
            % state and input constraints on tubes
            for ii = 1:T
                F = [F, alpha{ii} >= 0];
                for jj = 1:nJ
                   F = [F, ismember(X{ii}{jj}, Xc), ismember(U{ii}{jj}, Uc)];
                   if ii == T & ~isempty(Xf)
                      F = [F, ismember(X{ii}{jj}, Xf)]; 
                   end
                end
            end
            
            % constrain the states for {x_2, ..., x_T} 
            for ii_Delta_A  = 1:num_Delta_A
                for ii_Delta_B = 1:num_Delta_B
                   for ii = 1:T-1
                      for jj = 1:nJ
                        F = [F, Az*(Ahat + Delta_A_vertices{ii_Delta_A})*X{ii}{jj} + Az*(Bhat + Delta_B_vertices{ii_Delta_B})*U{ii}{jj}...
                             - Az*zc{ii+1} + max_W <= alpha{ii+1}*bz];
                      end          
                   end
                   % constrain the state for x_1
                   F = [F, Az*(Ahat + Delta_A_vertices{ii_Delta_A})*x0 + Az*(Bhat + Delta_B_vertices{ii_Delta_B})*u0 ...
                           - Az*zc{1} + max_W <= alpha{1}*bz];
                end
            end
            
%             % construct cost function
%             switch normType
%                 case 2
%                     cost_fcn = @(x,u) x'*Q*x + u'*R*u;
%                     cost_terminal = @(x) x'*P*x;
%                 case 1
%                     cost_fcn = @(x,u) Q(1)*norm(x, 1) + R(1)*norm(u, 1);
%                     cost_terminal = @(x) P(1)*norm(x, 1);
%                 case Inf
%                     cost_fcn = @(x,u) Q(1)*norm(x, Inf) + R(1)*norm(u, Inf);
%                     cost_terminal = @(x) P(1)*norm(x, Inf);
%                 otherwise
%                     warning('The norm type has to be 1, 2, or Inf.\n');
%             end
            
            cost_fcn = @(x,u) x'*Q*x + u'*R*u;
            cost_terminal = @(x) x'*P*x;
                    
            fprintf('Constructing objective function... \n');
            cost = cost_fcn(x0, u0);
            for ii = 1:T-1
                for jj = 1:nJ
                    cost = cost + cost_fcn(X{ii}{jj}, U{ii}{jj});
                end
            end
            % add terminal cost
            if ~isempty(P)
                for jj = 1:nJ
                   cost = cost + cost_terminal(X{T}{jj}); 
                end
            end
            cost = cost/(nJ*T);
            
            options = sdpsettings('solver','mosek', 'verbose', verbose);
            sol = struct;
            
            fprintf('Solver started...\n');
            
            if strcmp(type, 'value')
                diagnostics = optimize(F, cost);
            elseif strcmp(type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(F, cost, options, x0, U{1});
                return
            end
                        
            % extract the tubes
            tubes = cell(1, T);
            for ii = 1:T
                tubes{ii} = value(zc{ii}) + value(alpha{ii})*Z;
            end
            
            % extract the vertices
            vertices_set = cell(1, T);
            for ii = 1:T
               vertices = cell(1, nJ);
               for jj = 1:nJ
                  vertices{jj} = value(X{ii}{jj}); 
               end
               vertices_set{ii} = vertices;
            end
            
            % extract the control inputs
            u_tubes = cell(1,T);
            for ii = 1:T
                vertex_input = cell(1, nJ);
                for jj = 1:nJ
                    vertex_input{jj} = value(U{ii}{jj});
                end
               u_tubes{ii} = vertex_input; 
            end
            
            % save the solution
            sol.diagnostics = diagnostics;
            sol.solver_time = diagnostics.solvertime;
            sol.objective = value(cost);
            sol.u0 = value(u0);
            sol.x0 = x0;
            sol.tubes = tubes;   
            sol.u_tubes = u_tubes;
            sol.vertices_set = vertices_set;
            sol.status = diagnostics.problem;
            
            
            yalmip('clear');
        end

        %% Monimoy's method implementation
        function [sol] = SolveUniformDistFeedbackMPC(obj, type, verbose)
            % fixme: need to double check this
            
            if nargin < 2
                type = 'value';
                verbose = 2;
            elseif nargin < 3
                verbose = 2;
            end
            
            nx = obj.nx; nu = obj.nu;
            T = obj.horizon;
            
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            uncertain_system = obj.uncertain_system;
            
            Q = obj.Q; R = obj.R; Q_T = obj.Q_T;

           if strcmp(type, 'value')
                x0 = obj.x0;    
            elseif strcmp(type, 'optimizer')
                x0 = sdpvar(nx, 1);
            end
            
            % identify the augmented disturbance magnitude
            max_x_inf_norm = max(abs(obj.state_constr.V(:)));
            max_u_inf_norm = max(abs(obj.input_constr.V(:)));
            sigma_w_unif = eps_A*max_x_inf_norm + eps_B*max_u_inf_norm + sigma_w;
            
            % x = ZA x + ZB u + w           
            Z = kron(diag(ones(1,T),-1), eye(nx));
			ZA_block = Z*blkdiag(kron(eye(T), uncertain_system.A), zeros(nx, nx));
            ZB_block = Z*blkdiag(kron(eye(T), uncertain_system.B), zeros(nx, nu));
            
            % parameterize variables
            M = sdpvar( T*nu, T*nx, 'full');
            ubar = sdpvar(T*nu, 1);
            
            % updated
            M_aug = [zeros(T*nu, nx) M; zeros(nu, nx) zeros(nu, T*nx)];
            ubar_aug = [ubar; zeros(nu)];
            
            
            Id = eye(nx*(T+1));
            state_vec_robust = inv(Id- ZA_block)*(Id + ZB_block*M_aug);
            state_vec_robust_x0 = state_vec_robust(:, 1:nx);
            state_vec_robust_w = state_vec_robust(:, nx+1:end);
            state_vec_fixed = inv(Id- ZA_block)*(ZB_block*ubar_aug ) + state_vec_robust_x0*x0;
            
            % state vec = state_vec_fixed + state_vec_robust      
            
            constr = [];
            % structure constraint on M
            for ii = 1:T
               constr = [constr, M((ii-1)*nu+1:ii*nu, (ii-1)*nx+1:end) == zeros(nu, (T-ii+1)*nx)];
            end
            
            % robust state constraints
            Ax = obj.state_constr.A; bx = obj.state_constr.b;
            nFx = size(Ax, 1);
            try
                for ii = 1:T+1
                    for jj = 1:nFx
                    constr = [constr, Ax(jj, :)*state_vec_fixed((ii-1)*nx+1:ii*nx, 1) ...
                                            + norm(Ax(jj, :)*state_vec_robust_w((ii-1)*nx+1:ii*nx,:), 1)*sigma_w_unif <= bx(jj)];
                    end
                end
            catch ME
                sol = struct;
                sol.status = 1;
                sol.uniform_ub = sigma_w_unif;
                return
            end
            
            % robust terminal constraints
            if ~isempty(obj.terminal_constr)
                Af = obj.terminal_constr.A; bf = obj.terminal_constr.b;
                nFT = size(Af, 1);
                try
                    for jj = 1:nFT
                       constr = [constr, Af(jj, :)*state_vec_fixed(T*nx+1:(T+1)*nx, 1) ...
                                               + norm(Af(jj, :)*state_vec_robust_w(T*nx+1:(T+1)*nx, :), 1)*sigma_w_unif <= bf(jj)];
                    end
                catch
                    sol = struct;
                    sol.status = 1;
                    return
                end
            end
            % robust input constraint
            Au = obj.input_constr.A; bu = obj.input_constr.b;
            nFu = size(Au, 1);
            
            try
                for ii = 1:T
                    for jj = 1:nFu
                       constr = [constr, Au(jj,:)*ubar((ii-1)*nu+1:ii*nu,1) + norm(Au(jj,:)*M((ii-1)*nu+1:ii*nu, :), 1)*sigma_w_unif <= bu(jj)];
                    end
                end
            catch
                sol = struct;
                sol.status = 1;
                sol.uniform_ub = sigma_w_unif;
                return
            end
            % cost function
            cost = 0;
            Q_mat = kron(eye(T), Q);
            Q_mat = blkdiag(Q_mat, Q_T);
            cost = cost + state_vec_fixed'*Q_mat*state_vec_fixed;
            
            R_mat = kron(eye(T), R);
            cost = cost + ubar'*R_mat*ubar;
            
			ops = sdpsettings('verbose', verbose, 'solver', 'mosek');
      
            if strcmp(type, 'value')
                diagnostics = optimize(constr, cost, ops);
            elseif strcmp(type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(constr, cost, ops, x0, ubar(1:nu));
                return
            end
            
			if diagnostics.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(diagnostics.problem), '\n']);
                sol = struct;
                sol.diagnostics = diagnostics;
                sol.status = diagnostics.problem;
                sol.uniform_ub = sigma_w_unif;
                return
			end

			sol = struct;
            sol.M = value(M);
            sol.ubar = value(ubar);
            sol.cost = value(cost);
            sol.x0 = x0;
            sol.solver_time = diagnostics.solvertime;
            sol.diagnostics = diagnostics;
            sol.status = diagnostics.problem;
            sol.uniform_ub = sigma_w_unif;   
            yalmip('clear');
        end

        %% disturbance feedback
        function [sol] = SolveAugDistFeedbackSLSMPC(obj, type, opt)
            % the convex formulation for time-delay system does not exist:
            % we simply cannot convexify it. Instead, we shoud do use the
            % progressive approach to find the uncertainty norm bounds. 
            
            % we have constraints on u_T for this function
            
            if nargin < 2
                opt = struct;
                opt.verbose = 2;
                opt.solver = 'mosek';
                type = 'value';
            elseif nargin < 3
                opt = struct;
                opt.verbose = 2;
                opt.solver = 'mosek';
            end
            
            nx = obj.nx; nu = obj.nu; na = 0; nb = 0;
            T = obj.horizon;

            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            
            uncertain_system = obj.uncertain_system;
            
            if strcmp(type, 'value')
                x0 = obj.x0;    
            elseif strcmp(type, 'optimizer')
                x0 = sdpvar(nx, 1);
            end
            % construct system response variables 
			Phi_x = sdpvar( (T + 1) * nx, (T + 1) * nx, 'full');
			Phi_u = sdpvar( (T + 1) * nu, (T + 1) * nx, 'full');
            % extract matrix blocks
            Phi = [Phi_x; Phi_u];
  
            % construct the sigma matrix
            sigma_seq = sdpvar(1, T);
            Sigma_mat = eye(nx);
            for ii = 1:T
                Sigma_mat = blkdiag(Sigma_mat, sigma_seq(ii)*eye(nx));
            end
            
            % Construct the objective function using nominal cost
			sqrtQ = sqrtm(obj.Q);
			sqrtR = sqrtm(obj.R);
            sqrtQ_T = sqrtm(obj.Q_T);

			Qset = cell(T + 1, 1); Rset = cell(T + 1, 1);
            % no penalty on the initial state
			Qset{1} = zeros(size(sqrtQ));
			for i = 2 : T 
				Qset{i} = sqrtQ;
			end
			Qset{T + 1} = sqrtQ_T;

			for i = 1:T+1
				Rset{i} = sqrtR;
            end
            
            % attention: no penalty on u_T
			Rset{T + 1} = zeros(size(sqrtR));

			cost = 0;
            for i  = 1 : T+1
                % attention: add reference target state
                cost = cost + norm(Qset{i}*Phi_x((i - 1)*nx + 1: i*nx, 1:nx )*x0, 2)^2 ...
                       + norm(Rset{i}*Phi_u((i - 1)*nu + 1: i*nu, 1:nx )*x0, 2)^2;      
            end
            
            % add constraints to the problem
            constr = [];
            constr = [constr, sigma_seq >= 0];
            
            % Structure constraint of Phi_x and Phi_u
			for k = 1 : T
				constr = [constr, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (T + 1 - k)*nx)];
            end

			for k = 1 : T
				constr = [constr, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (T + 1 - k)*nx)];
            end

			Z = kron(diag(ones(1,T),-1), eye(nx));
			ZA_block = Z*blkdiag(kron(eye(T), uncertain_system.A), zeros(nx, nx));
            ZB_block = Z*blkdiag(kron(eye(T), uncertain_system.B), zeros(nx, nu));
			Id = eye((T + 1)*nx);
            % add affine constraint
            constr = [constr, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];

            % state constraints
			if ~isempty(obj.state_constr)
				Fx = obj.state_constr.A; bx = obj.state_constr.b;
				nFx = size(Fx, 1); nbx = length(bx);
            else
                warning('must have state constraints');
            end
            
            for ii = 1:T
               for jj = 1: nFx
                  f = Fx(jj,:); b = bx(jj);
                  LHS = f*Phi_x((ii-1)*nx+1:ii*nx,1:nx)*x0;
                  for kk = 1:ii-1
                      LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx),1);
                  end
                  constr = [constr, LHS <= b];  
               end
            end
            
            % terminal constraint   
            if ~isempty(obj.terminal_constr)
                Ft = obj.terminal_constr.A; bt = obj.terminal_constr.b;
                nFt = size(Ft, 1); nbt = length(bt);
            else 
                Ft = Fx; bt = bx; nFt = nFx; nbt = nbx;
            end
           
            
            for jj = 1:nFt
                f = Ft(jj,:); b = bt(jj);
                LHS = f*Phi_x(T*nx+1:(T+1)*nx,1:nx)*x0;
                for kk = 1:T
                   LHS = LHS + norm(f*Phi_x(T*nx+1:(T+1)*nx,kk*nx+1:(kk+1)*nx),1);
                end
                constr = [constr, LHS <= b];
            end
            
            % add input constraint
            if ~isempty(obj.input_constr)
				Fu = obj.input_constr.A; bu = obj.input_constr.b;
				nFu = size(Fu, 1); nbu = length(bu);
            else
                warning('must have input constraints');
            end
            
            for ii = 1:T+1
               for jj = 1: nFu
                  f = Fu(jj,:); b = bu(jj);
                  LHS = f*Phi_u((ii-1)*nu+1:ii*nu,1:nx)*x0;
                  for kk = 1:ii-1
                      LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),1);
                  end
                  constr = [constr, LHS <= b];  
               end
            end
            
            %  save inf_norm upper bounds for x_k
            x_upper_bounds = cell(1,T + 1);
            u_upper_bounds = cell(1,T + 1);

            for ii = 1:T+1
                x_ii_norm_bd = norm(Phi_x((ii-1)*nx+1:ii*nx, 1:nx)*x0, inf);
                for jj = 1:ii-1
                   x_ii_norm_bd = x_ii_norm_bd + norm(Phi_x((ii-1)*nx+1:ii*nx, jj*nx+1:(jj+1)*nx), inf);
                end
                x_upper_bounds{ii} = x_ii_norm_bd;
            end
            
            for ii = 1:T+1
               u_ii_norm_bd = norm(Phi_u((ii-1)*nu+1:ii*nu, 1:nx)*x0,inf); 
               for jj = 1:ii-1
                   u_ii_norm_bd = u_ii_norm_bd + norm(Phi_u((ii-1)*nu+1:ii*nu, jj*nx+1:(jj+1)*nx), inf);
               end
               u_upper_bounds{ii} = u_ii_norm_bd;
            end
            
            % constraints on sigma_0, sigma_1, ..., sigma_{T-1}
            for ii = 1:T
                x_upper_bd = x_upper_bounds{ii};
                u_upper_bd = u_upper_bounds{ii};
                constr = [constr, sigma_seq(ii) >= eps_A*x_upper_bd + eps_B*u_upper_bd + sigma_w];
            end
            
             % solve the problem
             if isfield(opt, 'solver')
                 solver = opt.solver;
             else
                 solver = 'mosek';
             end
             
             if isfield(opt, 'verbose')
                 verbose = opt.verbose;
             else
                 verbose = 2;
             end
             
			ops = sdpsettings('verbose', verbose, 'solver', solver);
            
            if strcmp(type, 'value')
                solution = optimize(constr, cost, ops);
            elseif strcmp(type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(constr, cost, ops, x0, Phi_u(1:nu, 1:nx));
                return
            end
                  
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
			end

			sol = struct;
			Phi_x_val = value(Phi_x); Phi_u_val = value(Phi_u);
            Phi_u_val(find(isnan(Phi_u_val)==1)) = 0;
            fval = value(cost);
            sigma_seq_val = value(sigma_seq);
            
            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            sol.fval = fval; 
            sol.sigma_seq = sigma_seq_val;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;

            obj.Phi_x = Phi_x_val; obj.Phi_u = Phi_u_val; obj.sigma_seq = sigma_seq_val;
            yalmip('clear');% clear all yalmip variables
        end
        
        %% simulate trajectories in prediction 
        function [traj_set] = SimulateTrajwithSystemResponses(obj, num_traj)
            % simulate the closed-loop trajectories in prediction using the
            % saved system responses Phi_x, Phi_u
            traj_set = cell(1, num_traj);
            
            uncertain_system = obj.uncertain_system;
            A_mat = uncertain_system.A;
            B_mat = uncertain_system.B;
            x0 = obj.x0;
            nx = obj.nx; nu = obj.nu; na = obj.na; nb = obj.nb;
            
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            
            T = obj.horizon;
            
            % synthesize the state feedback controller
            Phi_x = obj.Phi_x; Phi_u = obj.Phi_u;
            K = Phi_u*inv(Phi_x);
      
            for ii = 1:num_traj
               traj = struct;
               x_seq = [x0]; u_seq = []; w_seq = [];
               DeltaA_seq = cell(1, T); DeltaB_seq = cell(1, T);
               xsls = x0;
               
               Delta_A = DeltaOperator(nx, nx, 1, eps_A, 0);
               Delta_B = DeltaOperator(nx, nu, 1, eps_B, 0);
               for jj = 1:T
                   w = rand(nx,1).*2*sigma_w - sigma_w;           
                   w_seq = [w_seq w];

                   usls = zeros(nu, 1);
                   for kk = 1:jj
                       usls = usls + K((jj - 1)*nu + 1:jj*nu,(kk - 1)*nx + 1:kk*nx)*x_seq(:,kk);
                   end
                   u_seq = [u_seq usls];
                   
                   DeltaA_seq{jj} = Delta_A;
                   DeltaB_seq{jj} = Delta_B;
                   
                   xsls = (A_mat + Delta_A)*x_seq(:,end) + (B_mat + Delta_B)*u_seq(:,end) + w;
                   x_seq = [x_seq xsls];
                
               end
               
               traj.x_seq = x_seq;
               traj.u_seq = u_seq;
               traj.w_seq = w_seq;
               traj.DeltaA_seq = DeltaA_seq;
               traj.DeltaB_seq = DeltaB_seq; 
               traj_set{ii} = traj;
               
            end
                       
        end

        function [traj_set] = SimulateTrajwithTubeMPC(obj, tube_sol, num_traj)
            % recover model 
            T = obj.horizon; x0 = obj.x0;
            nx = obj.nx; nu = obj.nu; 
            
            uncertain_system = obj.uncertain_system;
            A = uncertain_system.A; B = uncertain_system.B;
 
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;
            u0 = tube_sol.u0;
            
            tube_seq = tube_sol.tubes;
            u_tubes = tube_sol.u_tubes;
            vertices_set = tube_sol.vertices_set;

            % extract Delta_A , Delta_B and w sequences from sls_traj_set
            DeltaA_set = cell(1,num_traj);
            DeltaB_set = cell(1,num_traj);
            w_set = cell(1, num_traj);

            % compute the tube MPC trajectories
            traj_set = cell(1, num_traj);
            for ii = 1:num_traj
               traj = struct;
               x_seq = [x0]; u_seq = [u0]; w_seq = [];
               DeltaA = DeltaOperator(nx, nx, 1, eps_A, 1);
               DeltaB = DeltaOperator(nx, nu, 1, eps_B, 1);
               DeltaA_set{ii} = DeltaA;
               DeltaB_set{ii} = DeltaB;
                
               xi_tube = x0; 
               for jj = 1:T
                   u_cur = u_seq(:, end);
                   A_robust = A + DeltaA;
                   B_robust=  B + DeltaB;

                   w = rand(nx,1).*2*sigma_w - sigma_w;           
                   w_seq = [w_seq w];

                   xi_tube = A_robust*xi_tube + B_robust*u_cur + w;

                   tube = vertices_set{jj};

                   [coef, status] = FindHullCoeff(xi_tube, cell2mat(tube)', 1);
                   if status ~= 0
                       keyboard;
                       traj_set = [];
                       warning('Tube MPC evaluation infeasible.');
                       return
                   end
                   
                   % find tube MPC control inputs
                   u_next = cell2mat(u_tubes{jj})*coef;
                   u_seq = [u_seq u_next];
                   x_seq = [x_seq xi_tube];

                   end
                   traj.xi_seq = x_seq;
                   traj.u_seq = u_seq(:, 1:end-1);
                   traj.w_seq = w_seq;
                   traj.DeltaA = DeltaA;
                   traj.DeltaB = DeltaB; 
                   traj_set{ii} = traj;
            end
 
        end
        
        function [traj_set] = SimulateTrajwithUniformDistFeedbackMPC(obj, df_sol, num_traj)
            % recover model 
            T = obj.horizon; x0 = df_sol.x0;
            nx = obj.nx; nu = obj.nu; 

            uncertain_system = obj.uncertain_system; 
            A_mat = uncertain_system.A; B_mat = uncertain_system.B;
            
            eps_A = obj.eps_A; eps_B = obj.eps_B; sigma_w = obj.sigma_w;

            % extract solution
            M = df_sol.M; ubar = df_sol.ubar;
            
            DeltaA_set = cell(1, num_traj); DeltaB_set = cell(1, num_traj);
            traj_set = cell(1, num_traj);
            for ii = 1:num_traj
               traj = struct;
               x_seq = [x0]; u_seq = []; w_seq = []; 
               xsls = x0;
               w_aug_seq = [];
               
               Delta_A = DeltaOperator(nx, nx, 1, eps_A, 1);
               Delta_B = DeltaOperator(nx, nu, 1 , eps_B, 1);
               DeltaA_set{ii} = Delta_A;
               DeltaB_set{ii} = Delta_B;
                   
               for jj = 1:T
                   w = rand(nx,1).*2*sigma_w - sigma_w;           
                   w_seq = [w_seq w];

                   usls = ubar((jj-1)*nu+1:jj*nu,:);
                   for kk = 1:jj - 1
                       usls = usls + M((jj - 1)*nu + 1:jj*nu,(kk - 1)*nx + 1:kk*nx)*w_aug_seq(:,kk);
                   end
                   u_seq = [u_seq usls];
                   xsls = (A_mat + Delta_A)*x_seq(:,end) + (B_mat + Delta_B)*u_seq(:,end) + w;
                   xsls_nominal = A_mat*x_seq(:,end) + B_mat*u_seq(:,end);
                   w_aug = xsls - xsls_nominal;

                   x_seq = [x_seq xsls];
                   w_aug_seq = [w_aug_seq w_aug];
               end

               traj.x_seq = x_seq;
               traj.u_seq = u_seq;
               traj.w_seq = w_seq;
               traj.w_aug_seq = w_aug_seq;
               traj.Delta_A = Delta_A;
               traj.Delta_B = Delta_B; 
               traj_set{ii} = traj;

            end
        end
        
    end
    
end

