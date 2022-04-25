classdef UncertainLTISystem < handle
    %UNCERTAINLTISYSTEM class of an LTI system with uncertainty
    %   x(k+1) = Ax(k) + Bu(k) + D_A x(k) + D_B u(k) + w(k)
    
    properties(SetAccess = public)
        A; B; % system dynamics matrices
        nx; nu; % state and input dimension
        x0; % initial state
        sigma_w = 0; % the ell_inf norm of additive disturbances w(k)
        eps_A = 0; eps_B = 0; % ||D_A||_inf <= eps_A, ||D_B||_inf <= eps_B
        DA_vertices; DB_vertices; % cell of vertices of the model uncertainty
        horizon = 0; % horizon of a optimal control problem, used to generate relevant block matrices
        K = []; % u = Kx is the local linear feedback controller
        OLRIS; CLRIS; % open-loop and closed-loop (with local controller u = Kx) robust invariant set
    end
    
    methods
        function obj = UncertainLTISystem(params)
            %UNCERTAINLTISYSTEM Construct an instance of this class
            %   initialize with a given struct params
 
            obj.A = params.A;
            obj.B = params.B;
            obj.x0 = params.x0;
            
            nx = size(params.A, 1);
            nu = size(params.B, 2);
            obj.nx = nx;
            obj.nu = nu;
            
            if isfield(params, 'sigma_w')
                obj.sigma_w = params.sigma_w;
            end

            if isfield(params, 'eps_A')
                obj.eps_A = params.eps_A;
            end
            
            if isfield(params, 'eps_B')
                obj.eps_B = params.eps_B;
            end
            
            if isfield(params, 'horizon')
                obj.horizon = params.horizon;
            end
            
            D_A_dim = size(obj.A); D_B_dim = size(obj.B);
            
            if isfield(params, 'Delta_A_vertices')
                obj.DA_vertices = params.Delta_A_vertices;
            else
                obj.DA_vertices = FindMatVertices(D_A_dim, obj.eps_A);
            end
            
            if isfield(params, 'Delta_B_vertices')
                obj.DB_vertices = params.Delta_B_vertices;
            else
                obj.DB_vertices = FindMatVertices(D_B_dim, obj.eps_B);
            end
           
        end
        
        %% find local controller u = Kx through LQR
        function [K, P] = find_K_LQR(obj, Q, R)
            [K, P] = dlqr(obj.A, obj.B, Q, R);
            obj.K = -K;
        end
        
        %% find minimum disturbance invariant set Z_inv
        function [output, isConverge] = minInvSet(obj, N)
            % N: maximum number of iterations
            sigma_w = obj.sigma_w;
            W0 = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));;
            A = obj.A;
            B = obj.B;
            n = size(A, 1); m = size(B, 2);
            K = obj.K;
            if isempty(K)
                error('Local controller u = Kx has to be computed first.');
            end
            Acl = A + B*K;
            W = W0;
            isConverge = 0;
            volList = zeros(1, N+1);
            volList(1) = W.volume();
            for ii = 1:N
               W = W + Acl^ii*W0; % find disturbance invariant set using iterations
               W.minHRep(); % Important: call minHRep to reduce redundant inequalities
               volList(ii+1) = W.volume();
               if abs(volList(ii+1) - volList(ii)) <1e-2
                   isConverge = 1;
                   break;
               end
            end
            output = W;
        end
        
        %% preset
        function [PreS] = preset(obj, Xc, Uc)
           % compute the preset of Xc for all possible u \in Uc, using
           % projection methods
          
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;

            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 

            Fbar = [F*A F*B; zeros(nU, n) G]; fbar = [f; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);

            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims);
        end

        %% presetRobust
        function [PreS] = presetRobust(obj, Xc, Uc)
            % compute the preset of Xc for all possible u \in Uc
            W = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            robustXc = Xc - W;
            F = robustXc.A; f = robustXc.b;
            
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;

            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 
               
            numDA_vertices = length(DA_vertices); numDB_vertices = length(DB_vertices);
            
            Fbar = []; fbar = [];
            
            % in case of no model uncertainty
            constr = [F*A F*B];
            Fbar = [Fbar; constr]; fbar = [fbar; f];
                    
            for ii = 1:numDA_vertices
                for jj = 1:numDB_vertices
                    constr = [F*(A + DA_vertices{ii}) F*(B + DB_vertices{jj})];
                    Fbar = [Fbar; constr]; fbar = [fbar; f];
                end
            end
                
            newRow = [zeros(nU, n) G];
            Fbar = [Fbar; newRow]; fbar = [fbar; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);
            tic
            fprintf('remove lift poly redundancy.\n');
            liftedPolyhedron.minHRep();      
            toc
            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims); 
            tic
            fprintf('remove preS redundancy.\n');
            PreS.minHRep();
            toc
        end
        
        %% preset of autonomous system
        function [PreS] = presetAuto(obj, Xc, sysA)
            % find preset of autonomous system x_+  = sysA*x
            F = Xc.A;  f= Xc.b;
            PreS = Polyhedron(F*sysA, f);
        end
        
        %% preset of autonomous system under model uncertainty
        function [PreS] = presetAutoRobust(obj, Xc)
            % Enumeration of vertices of the uncertainty set is applied.
            W = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            robustXc = Xc - W;
            F = robustXc.A; f = robustXc.b;

            A = obj.A; B = obj.B; 
            DA_vertices = obj.DA_vertices; DB_vertices = obj.DB_vertices;
            
            % use the local controller
            K = obj.K;
            
            n = size(A,2); m = size(B, 2);
            numDA_vertices = length(DA_vertices); numDB_vertices = length(DB_vertices);
            
            Fbar = []; fbar = [];
            for ii = 1:numDA_vertices
                for jj = 1:numDB_vertices
                    % fixme: A+BK is used since K = obj.K
                    constr = [F*((A + DA_vertices{ii})+(B + DB_vertices{jj})*K)];
                    Fbar = [Fbar; constr]; fbar = [fbar; f];
                end
            end
            
            PreS = Polyhedron(Fbar, fbar);
 
        end

        %% robustInvariantset
        function [RIS, diagnostic] = robustInvariantSet(obj, Xinit, Uc, Nstep, options)
            % Xc is the given initial set on states; Uc is the constraint
            % on u; Nstep is the maximum step simulated forward.
            % This function computes the robust control invariant set.
            if nargin < 4
                Nstep = 10;
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            elseif nargin < 5
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            end
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            volList = zeros(1, Nstep+1);
            volList(1) = Xinit.volume();
            for ii = 1:Nstep
                fprintf('Invariant set iter %d/%d \n', ii, Nstep);
                fprintf('volume: %f \n', volList(ii));
                
                diagnostic.runningStep = ii;
               
%                 figure; Xinit.plot(); 
%                 title(['step = ', num2str(ii), ' total step = ', num2str(Nstep)]);
%                 pause(0.5);
                if options.robust == 0
                    preS = obj.preset(Xinit, Uc);
                else 
                    preS = obj.presetRobust(Xinit, Uc);
                end
                X_new = and(Xinit, preS);
                X_new.minHRep();
                
                save('RIS_temp', 'X_new');
                
                if X_new == Xinit
                    diagnostic.converge = 1;
                    break;
                end
                Xinit = X_new;
                volList(ii+1) = Xinit.volume();
                
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
            end
            RIS = Xinit;
            obj.OLRIS = RIS;
            diagnostic.volList = volList;
        end
        
        %% Robust invariant set for closed-loop systems
        function [RIS, diagnostic] = robustInvariantSetClosedLoop(obj, Xc, Uc, Nstep, options)
            % compute the robust invariant set of the closed-loop system
            % A+BK
            if nargin < 5
               options = struct;
               options.plot = 0;
               options.robust = 0;
               options.minVol = 0.5;
            end
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;
            
            A = obj.A; B = obj.B; 
            K = obj.K;
            
%             [F_unify, G_unify, ~] = convert_Poly2Mat(Xc, Uc);
%             Xc_new = Polyhedron(F_unify + G_unify*K, ones(size(F_unify, 1), 1));
            Fbar = [F; G*K]; fbar = [f; g];
            Xc_new = Polyhedron(Fbar, fbar);
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            Xinit = Xc_new;
            sysA = A + B*K;
            
            volList = zeros(1, Nstep+1);
            volList(1) = Xinit.volume();
            for ii = 1:Nstep
                fprintf('Invariant set iter %d/%d \n', ii, Nstep);
                diagnostic.runningStep = ii;
                if options.plot ~= 0
                    clf;
                    Graphics.show_convex(Xinit, 'r');
                    if isfield(options, 'xrange') && ~isempty(options.xrange)
                         xlim(options.xrange);
                    end
                    if isfield(options, 'yrange') && ~isempty(options.yrange)
                        ylim(options.yrange);
                    end
%                     Xinit.plot(); 
                    title(['step = ', num2str(ii), ' total step = ', num2str(Nstep)]);
                    pause(0.5);
                end
                if options.robust ~= 1
                    PreS = obj.presetAuto(Xinit, sysA);
                else 
                    PreS = obj.presetAutoRobust(Xinit);
                end
                X_new = and(Xinit, PreS);
                volList(ii+1) = X_new.volume();
                if X_new == Xinit
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                Xinit = X_new;
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
            end
            RIS = Xinit;
            obj.CLRIS = RIS;
            diagnostic.volList = volList;
        end
        
        %% compute local stabilizing controller u = K xi
        function [K_value] = local_stabilizing_controller(obj)
           nx = obj.nx; nu = obj.nu;           
           A = obj.A;
           B = obj.B;
           
           G = sdpvar(nx, nx);
           L = sdpvar(nu, nx);
           
           rho = 0.99;
           constr = [];
           
           DA_vertices = obj.DA_vertices;
           DB_vertices = obj.DB_vertices;

            num_A_vert = length(DA_vertices);
            num_B_vert = length(DB_vertices);
            
            for ii = progress(1:num_A_vert)
                for jj = 1:num_B_vert
                    A_robust = A + DA_vertices{ii};
                    B_robust = B + DB_vertices{jj};
                    mat = [rho*G (A_robust*G + B_robust*L)';
                           A_robust*G+B_robust*L  G];
                    constr = [constr, mat >= 0];
                end
            end
           
           optimize(constr);
           G_value = value(G);
           L_value = value(L);
           K_value = L_value*inv(G_value);
           obj.K = K_value;
            
        end
        
    end
end

