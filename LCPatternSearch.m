function [x_hist, f_hist] = LCPatternSearch(x0, obj_fun, A, l, u, sl_ini, cvg_par, se_par)
    %==========================================================================
    %Yang, JS; 2020-08-09
    %A pattern search algorithm for linearly constrained optimization
    %==========================================================================
    %INPUT LIST
    %--------------------------
    %   x0     : initial point
    %   obj_fun: function handle of objective
    %   A, l, u: define linear constraints l <= Ax <= u
    %   sl_ini : initial step length parameter of pattern
    %   cvg_par: a vector to define convergent conditions
    %       sl_esp  = cvg_par(1): the convergent parameter for step length of pattern
    %       ite_max = cvg_par(2): the maximum iteration number for pattern search
    %   se_par : parameters for the adjustion of the step length parameter
    %       shrink = 0< se_par(1) < 1
    %       extend =    se_par(2) > 1
    
    %OUTPUT LIST
    %--------------------------
    %   x_hist : the sequence of intermediate points
    %   f_hist : the sequence of objective values at intermediate points
    %==========================================================================

    %number of variables
    n_x     = length(x0); 

    %convergent parameters
    sl_esp  = cvg_par(1);
    ite_max = cvg_par(2);

    %initialization
    sl      = sl_ini;
    ite     = 0;
    x       = x0;
    f       = obj_fun(x);
    x_hist  = zeros(n_x, ite_max);
    f_hist  = zeros(1, ite_max);
    fprintf('------------------------------------------------------------->\n')
    fprintf('Initial design:\n');
    valuedisplay(x, 'x', 5)
    fprintf('------------------------------------------------------------->\n')
    fprintf('Initial function value:\n');
    fprintf('f_val = %15.6f\n', f);
    fprintf('------------------------------------------------------------->\n')

    x_hist(:,1) = x;
    f_hist(:,1) = f;

    %start iteration
    while sl >= sl_esp && ite <= ite_max

        ite = ite + 1;

        fprintf('=============================================================>\n')
        fprintf('Iteration number: %d\n', ite);
        fprintf('=============================================================>\n')

        eps0 = 1e-3;
        patternMat = pattern(A, l, u, x, eps0);

        [x, f, flag] = polls(x, f, sl, A, l, u, patternMat, obj_fun);

        %record the point and function value
        x_hist(:,ite+1) = x;
        f_hist(:,ite+1) = f;
        fprintf('------------------------------------------------------------->\n')
        fprintf('Current design:\n');
        valuedisplay(x, 'x', 5);
        fprintf('------------------------------------------------------------->\n')
        fprintf('Current function value:\n');
        fprintf('f_val = %15.6f\n', f);
        fprintf('------------------------------------------------------------->\n')

        %update the parameters of step length
        sl = slUpdate(sl, flag, se_par);
    end

    fprintf('------------------------------------------------------------->\n')
    fprintf('The step length parameter is %15.6f\n', sl);
    fprintf('=============================================================>\n')
    if ite >= ite_max
        fprintf('Pattern Search terminated since the iteration number reached the maximum iteration number!\n');
    else
        fprintf('Pattern Search terminated since the step length parameter is less than the tolerence!\n')
    end
    fprintf('=============================================================<\n')

    x_hist = x_hist(:,1:ite+1);
    f_hist = f_hist(:,1:ite+1);

end

%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: define the pattern
%-----------------------------------------------
function patternMat = pattern(A, l, u, xk, eps0)
    %Yang, JS; 2020-08-09

    %generator set
    %-------------------------
    %number of linear constraints
    n_lc = length(l);
    %number of variable
    n_x  = length(xk);

    I_l  = zeros(1, n_lc);
    I_u  = zeros(1, n_lc);

    fullColumnRank = 0;
    i_fcrk = 0;

    dist_l = zeros(1,n_lc);
    dist_u = zeros(1,n_lc);
    for ii = 1:1:n_lc
        %distance between xk and the boundary of the feasible domain
        dist_l(ii) = abs(A(ii,:)*xk-l(ii)) / norm(A(ii,:));
        dist_u(ii) = abs(A(ii,:)*xk-u(ii)) / norm(A(ii,:));
    end

    while ~fullColumnRank
        i_fcrk = i_fcrk+1;
        epsk = 0.5^(i_fcrk-1)*eps0;
        I_u(dist_u <= epsk) = 1;
        I_l(dist_l <= epsk) = 1;

        v_l = -A(I_l == 1, :)';
        u_l =  A(I_u == 1, :)';
        V = [v_l, u_l];

        if rank(V) == sum(I_u)+sum(I_l)
            fullColumnRank = 1;
        end
    end

    %rational positive basis
    %-------------------------
    V(:,(~any(V))) = [];
    VV   = V'*V;
    II0 = eye(size(VV));
    II1 = eye(n_x);
    VVV  = V*(VV\II0);
    VVVV = VVV*V';

    % N    = VVVV;
    N     = rref((II1-VVVV)')';
    N     = N(:, sum(abs(N))>0);

    %Pattern Matrix
    N2    = [N, -N];
    VVV2  = [VVV, -VVV];
    patternMat = [N2, VVV2];

end
%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: poll for pattern directions
%-----------------------------------------------
function [x, f, flag] = polls(x0, f0, sl, A, l, u, patternMat, obj_fun)
    %Yang, JS; 2020-08-08

    %poll flag: 1 for successful poll; 0 for unsuccessful poll
    flag       = 0;
    %number of variables
%     n_x        = length(x0);

    %number of pattern direction
    n_d        = size(patternMat,2);
    %trial and polls
    x = x0;
    f = f0;
    f_trial = f0;
    ctr_trial = 0;
    for ii = 1:1:n_d
        %trial points
        x_trial   = x0 + sl*patternMat(:,ii);
        if sum(A*x_trial > u) == 0 && sum(A*x_trial < l) ==0
            %function value at trial points
            f_trial   = obj_fun(x_trial);
            ctr_trial = ctr_trial+1;
        end
        
        if f_trial < f0
            x = x_trial;
            f = f_trial;
            flag = 1;
            ctr_trial = ctr_trial+1;
            break;
        end

    end

    %output the number of trial
    fprintf('------------------------------------------------------------->\n')
    fprintf('Number of trial in the current step is %d\n', ctr_trial);
    fprintf('------------------------------------------------------------->\n')
    
end
%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: update of step length parameter
%-----------------------------------------------
function sl = slUpdate(sl, flag, se_par)
    %Yang, JS; 2020-08-08

    shrink = se_par(1);
    extend = se_par(2);

    if flag == 0
        sl = sl*shrink;
    else
        sl = sl*extend;
    end

end
%==========================================================================