function [x_hist, f_hist] = UCPatternSearch(x0, obj_fun, sl_ini, cvg_par, se_par, pattern)
    %==========================================================================
    %Yang, JS; 2020-08-08
    %A pattern search algorithm for unconstrained optimization
    %==========================================================================
    %INPUT LIST
    %--------------------------
    %   x0     : initial point
    %   obj_fun: function handle of objective
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
    fprintf('Initial design:\n');
    valuedisplay(x, 'x', 5)
    fprintf('Initial function value:\n');
    fprintf('f_val = %15.6f\n', f);

    x_hist(:,1) = x;
    f_hist(:,1) = f;

    %start iteration
    while sl >= sl_esp && ite <= ite_max

        ite = ite + 1;

        fprintf('=============================================================>\n')
        fprintf('Iteration number: %d\n', ite);
        fprintf('=============================================================>\n')

        if strcmpi(pattern, 'HJ')
            if ite == 1
                [x, f, flag] = polls(x, f, sl, obj_fun);
            else
                zooms  = x_hist(:,ite) - x_hist(:,ite-1);
                x_tmpc = x+zooms;
                x_tmp0 = x;
                [x, f, flag] = polls(x_tmpc, f, sl, obj_fun);
                if flag == 0
                    [x, f, flag] = polls(x_tmp0, f, sl, obj_fun);
                end
            end

        elseif strcmpi(pattern, 'CS')
            [x, f, flag] = polls(x, f, sl, obj_fun);
        else
            error('The pattern can only be <CS>(Coordinate Search) or <HJ>(Hooke and Jeeves)!');
        end

        %record the point and function value
        x_hist(:,ite+1) = x;
        f_hist(:,ite+1) = f;
        fprintf('Current design:\n');
        valuedisplay(x, 'x', 5);
        fprintf('Current function value:\n');
        fprintf('f_val = %15.6f\n', f);

        %update the parameters of step length
        sl = slUpdate(sl, flag, se_par);
    end

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
%SUBFUNCTION FOR PATTERN SEARCH: poll for pattern directions
%-----------------------------------------------
function [x, f, flag] = polls(x0, f0, sl, obj_fun)
    %Yang, JS; 2020-08-08

    %poll flag: 1 for successful poll; 0 for unsuccessful poll
    flag       = 0;
    %number of variables
    n_x        = length(x0);
    %matrix of pattern direction (positive,half)
    patternMat = eye(n_x);

    %trial and polls
    x = x0;
    f = f0;
    f_trial_min = f0;
    ctr_trial = 0;
    for ii = 1:1:n_x
        %trial points
        x_trial   = x + sl*patternMat(:,ii);
        %function value at trial points
        f_trial   = obj_fun(x_trial);
        f_trial_min = min(f_trial_min, f_trial);
        ctr_trial = ctr_trial + 1;

        if f_trial >= f
            %trial points
            x_trial   = x - sl*patternMat(ii);
            %function value at trial points
            f_trial   = obj_fun(x_trial);
            f_trial_min = min(f_trial_min, f_trial);
            ctr_trial = ctr_trial + 1;
        end
        
        if f_trial < f
            x = x_trial;
            f = f_trial;
        end

    end

    if f_trial_min < f0
        flag = 1;
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