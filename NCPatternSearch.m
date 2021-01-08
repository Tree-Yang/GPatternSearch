function [x_hist, f_hist, c_hist] = NCPatternSearch(x0, obj_fun, nlcon_fun, bounds, A, l, u,...
                                                sl_ini, cvg_par, se_par, pattern)
    %==========================================================================
    %Yang, JS; 2020-08-10
    %A pattern search algorithm for nonlinearly constrained optimization
    %==========================================================================
    %INPUT LIST
    %--------------------------
    %   x0     : initial point
    %   obj_fun: function handle of objective: min obj_fun(x)
    %   nlcon_fun: nonlinear constraints (vector): nlcon_fun(x) <= 0
    %   bounds : upper and lower bounds of variables: lb <= x <= ub
    %       lb = bounds(:,1): the lower bounds of variables
    %       ub = bounds(:,2): the upper bounds of variables
    %   A, l, u: Matrix and vectors to define linear constraint: l <= Ax <=u 
    %   sl_ini : initial step length parameter of pattern
    %   cvg_par: a vector to define convergent conditions
    %       sl_esp  = cvg_par(1): the convergent parameter for step length of pattern
    %       ite_max = cvg_par(2): the maximum iteration number for pattern search
    %   se_par : parameters for the adjustion of the step length parameter
    %       shrink = 0< se_par(1) < 1
    %       extend =    se_par(2) > 1
    %   pattern: type of pattern for pattern search
    %       CS : coordinate search
    %       HJ : pattern of Hooke and Jeeves
    
    %OUTPUT LIST
    %--------------------------
    %   x_hist : the sequence of intermediate points
    %   f_hist : the sequence of objective values at intermediate points
    %   c_hist:  the sequence of nonlinear constraint at intermediate points
    %==========================================================================

    %==========================================================================
    %ALGORITHM SETTING
    %-----------------------------------------------------
    %parameters for bisection procedure
    bisecn_max  = 10;               %maximum number of bi-section step
    bisec_eps   = 1e-3;             %tolorence of bi-section step
    %parameters of convergence condition
    delta_eps   = 1e-3;             %tolorence of step length of pattern search
    eta_eps     = 1e-3;             %tolorence of nonlinear constraints
    ite_max     = 100.0;            %maximum number of iterations
    global nTrialGPS                %a global counter of function evaluations
    nTrialGPS   = 0;
    
    %INITIALIZATION
    %-----------------------------------------------------
    %initial value of nonlinear constraints
    nlc0      = nlcon_fun(x0);
    %initial value of objective function
    f0        = obj_fun(x0);
    %number of variables
    n_x       = length(x0);
    %number of nonlinear constraints
    n_nlc     = length(nlc0);
    nTrialGPS = nTrialGPS + 1;
    
    %update parameters of parameters
    %The following parameters are defined accoording to the paper named "A global convergent
    %   augmented Lagrangian algorithm for optimization with general constraints and
    %   simple bounds", by Conn AR, Gould NI and Toint PL, 1991.
    eta0        = 1.0;
    omega0      = 1.0;
    alpha_eta   = 0.1;
    alpha_omega = 1.0;
    beta_eta    = 0.9;
    beta_omega  = 1.0;
%     gamma1  = min((0.1/omega0)^(1/alpha_omega), (0.1/eta0)^(1/alpha_eta));
%     gamma1  = min(1.0, gamma1);
    gamma1      = 0.1;

    %initial penalty parameters
    mu0         = 10.0;
    tau         = 0.01;
    mu_k        = mu0;

    %scaling matrix
    S           = eye(n_nlc);
    
    %initial estimate of multipliers (copy from matlab built-in function)
    shiftConst  = 1e-3;
    shift       = max(0,nlc0) + shiftConst;
    alphaL      = min(1,alpha_eta/(1-alpha_eta));
    lambda0     = (shift./mu0).^(1/alphaL);
    lambda0     = max(zeros(n_nlc,1), lambda0+diag(S).*nlc0./mu_k);
    lambda_k    = lambda0;
    
    alpha_k     = min(mu_k, gamma1);
%     omega_k     = min(0.1,omega0 * alpha_k^alpha_omega);
    omega_k     = omega0 * alpha_k^alpha_omega;
%     meshConst   = 1e3;
    theta_k     = 1.0/(1.0+norm(lambda0)+1.0/mu0);
%     delta_k     = meshConst * theta_k * omega_k;
    delta_k     = theta_k * omega_k;
    delta_k     = min(delta_k, sl_ini);
%     eta_k       = min(0.1,eta0 * alpha_k^alpha_eta);
    eta_k       = eta0 * alpha_k^alpha_eta;
    
    %initialize the decision variables, constraint function value, and nonlinear constraint values    
    x           = x0;
    f           = f0;
    nlc         = nlc0;
    
    %initialze feasible descision variables for bi-section procedure
    x_old       = x0;
    f_old       = f0;
    nlc_old     = nlc0;
    
    %initialize the stack of iteration history
    x_hist      = zeros(n_x, ite_max+1);
    f_hist      = zeros(1, ite_max+1);
    c_hist      = zeros(n_nlc, ite_max+1);
    x_hist(:,1) = x0;
    f_hist(:,1) = f0;
    c_hist(:,1) = nlc0;

    startTime    = datestr(now);
    fprintf('***************************************************************************************\n');
    fprintf('\t\t\t\t\t\t\t Algorithm Initialized!\n');
    fprintf('Start at: %s\n', startTime);
    %print initial points and function values
    fprintf('Initial point:\n');
    valuedisplay(x0, 'x', 5);
    fprintf('Objective of initial point:\n');
    fprintf('f_o = %15.6f\n', f0);
    fprintf('Nonlinear constraints of initial point:\n');
    valuedisplay(nlc0, 'nlc', 5);
    fprintf('***************************************************************************************\n');
    tic;

    %ITERATION OF SUB-PROBLEMS
    %-----------------------------------------------------
    for ite = 1:1:ite_max

        fprintf('===================================================================================>\n');
        fprintf('\t\t\t\t Iteration Number:\t %d\n', ite);
        fprintf('===================================================================================>\n');

        if isempty(A)
        %solve the bound-constrained sub problem
            [x_hist0, f_hist0, nlc_hist0, alf_hist0] = BCPatternSearchSub(x, f, nlc, obj_fun, ...
                                    bounds, nlcon_fun, sl_ini, S, lambda_k, mu_k ,cvg_par,...
                                    se_par, pattern);
        else
            %solve the bound-constrained sub problem
            [x_hist0, f_hist0, nlc_hist0, alf_hist0] = LCPatternSearchSub(x, f, nlc, obj_fun, ...
                                    A, l, u, nlcon_fun, sl_ini, S, lambda_k, mu_k ,cvg_par,...
                                    se_par);
        end

        %save the iteration history of the sub-problem to files
        subproblem_resfn = ['SubOptimization_',num2str(ite),'.mat'];
        save(subproblem_resfn, 'x_hist0', 'f_hist0', 'nlc_hist0', 'alf_hist0', '-v7.3');

        %result of the sub-problem
        x_sub = x_hist0(:, end);
        f_sub = f_hist0(:, end);
        nlc_sub = nlc_hist0(:, end);
        %print solution of every sub-problems
        fprintf('Solution of sub-problem-%d:\n', ite);
        valuedisplay(x_sub, 'x', 5);
        fprintf('Objective of sub-problem-%d:\n', ite);
        fprintf('f_o = %15.6f\n', f_sub);
        fprintf('Nonlinear constraints of sub-problem-%d:\n', ite);
        valuedisplay(nlc_sub, 'nlc', 5);

        %feasibility adaption
        %+++++++++++++++++++++++++++++
        if sum(nlc_sub>0) ~= 0
            fprintf('===================================================================================>\n');
            fprintf(['The solution of the sub-problem is not feasible. ',...
                        'A bisection strategy is activated!\n']);
            fprintf('===================================================================================>\n');
            %bisection
            %+++++++++++++++++++++++++++++
            [x_fea,f_fea,nlc_fea] = bisec(ite, x_sub, f_sub, nlc_sub, x_old, f_old, nlc_old, obj_fun, ...
                                            nlcon_fun, bisecn_max, bisec_eps);
            x = x_fea;
            f = f_fea;
            nlc = nlc_fea;

        end

        %update and record the variables, objective and constraint functions
        %+++++++++++++++++++++++++++++
        %feasible decision variables for bi-section
        x_old = x;
        f_old = f;
        nlc_old = nlc;
        %initial decision variables for next iteration
        x = x_sub;
        f = f_sub;
        nlc = nlc_sub;
        x_hist(:, ite+1) = x;
        f_hist(:, ite+1) = f;
        c_hist(:, ite+1) = nlc;

        if norm(nlc) <= eta_k
            %convergence check
            %+++++++++++++++++++++++++++++
            if delta_k <= delta_eps && norm(nlc) <= eta_eps
                fprintf('===================================================================================>\n');
                fprintf(['The algorithm converged since the mesh size ',...
                            'of the pattern is small enough!\n']);
                fprintf('===================================================================================>\n');
                break;
            else
                %parameters update (Lagrangian multiplier)
                %+++++++++++++++++++++++++++++
                lambda_k = max(zeros(n_nlc,1), lambda_k+diag(S).*nlc./mu_k);
                alpha_k  = min(mu_k, gamma1);
%                 omega_k  = min(0.1, omega_k * alpha_k^beta_omega);
                omega_k  = omega_k * alpha_k^beta_omega;
                theta_k  = 1/(1+norm(lambda_k)+1/mu_k);
                delta_k  = omega_k * theta_k;
%                 eta_k    = min(0.1, eta_k * alpha_k^beta_eta);
                eta_k    = eta_k * alpha_k^beta_eta;
            end
        else
            %parameters update (panalty parameters)
            %+++++++++++++++++++++++++++++
            mu_k     = tau * mu_k;
            alpha_k  = min(mu_k, gamma1);
%             omega_k  = min(0.1, omega0 * alpha_k^alpha_omega);
            omega_k  = omega0 * alpha_k^alpha_omega;
            theta_k  = 1/(1+norm(lambda_k)+1/mu_k);
            delta_k  = omega_k * theta_k;
%             eta_k    = min(0.1, eta0 * alpha_k^alpha_eta);
            eta_k    = eta0 * alpha_k^alpha_eta;
        end

    end

    if ite == ite_max
        fprintf('===================================================================================>\n');
        fprintf(['The algorithm converged since the iteration ',...
                    'number reached the maximum iteration number!\n']);
        fprintf('===================================================================================>\n');
    end

    x_hist      = x_hist(:,1:ite+1);
    f_hist      = f_hist(:,1:ite+1);
    c_hist      = c_hist(:,1:ite+1);
    
    save('NCPS_hist.mat', 'x_hist', 'f_hist', 'c_hist', '-v7.3');
    
    endTime    = datestr(now);
    ts = toc;

    fprintf('***************************************************************************************\n');
    fprintf('\t\t\t\t\t\t\t Algorithm Converged!\n');
    fprintf('Terminate at: %s\n', endTime);
    fprintf('Total number of function evaluation is: %d.\n',nTrialGPS);
    fprintf('Time consumption is: %15.6f s.\n', ts);
    %print final points and function values
    fprintf('Initial point:\n');
    valuedisplay(x, 'x', 5);
    fprintf('Objective of initial point:\n');
    fprintf('f_o = %15.6f\n', f);
    fprintf('Nonlinear constraints of initial point:\n');
    valuedisplay(nlc, 'nlc', 5);
    fprintf('***************************************************************************************\n');

end