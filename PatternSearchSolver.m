function [x_hist, f_hist, c_hist] = PatternSearchSolver(x0, obj_fun, bounds, ...
                                                            A, l, u, nlcon_fun,options)
    
    %Yang, Jiashu; 2020-08-11

    %x0       : initial point
    %obj_fun  : objective function
    %bounds   : bound constraints
    %A, l, u  : linear constraints
    %nlcon_fun: nonlinear constraints
    %options  : other parameters of the algorithm

    c_hist = [];

    %select proper solver for different input
    if isempty(bounds) && isempty(A) && isempty(l) ...
        && isempty(u) && isempty(nlcon_fun)
        %unconsitrained pattern search
        [x_hist, f_hist] = UCPatternSearch(x0, obj_fun, options.sl_ini,...
                            options.cvg_par, options.se_par, options.pattern);
    elseif ~isempty(bounds) && isempty(A) && isempty(l) ...
            && isempty(u) && isempty(nlcon_fun)
        %bound constrained pattern search
        [x_hist, f_hist] = BCPatternSearch(x0, obj_fun, bounds, options.sl_ini,...
                            options.cvg_par, options.se_par, options.pattern);
    elseif isempty(bounds) && ~isempty(A) && ~isempty(l) ...
            && ~isempty(u) && isempty(nlcon_fun)
        %linearly constrained pattern search
        [x_hist, f_hist] = LCPatternSearch(x0, obj_fun, A, l, u,...
                            options.sl_ini, options.cvg_par, options.se_par);
    elseif ~isempty(bounds) && isempty(A) && isempty(l) ...
            && isempty(u) && ~isempty(nlcon_fun)
        %nonlinearly and bound constrained pattern search
        [x_hist, f_hist, c_hist] = NCPatternSearch(x0, obj_fun, nlcon_fun, bounds,[], [], [],...
                            options.sl_ini, options.cvg_par, options.se_par, options.pattern);
    elseif  isempty(bounds) && ~isempty(A) && ~isempty(l) ...
        && ~isempty(u) && ~isempty(nlcon_fun)
        %nonlinearly and linearly constrained pattern search
        [x_hist, f_hist, c_hist] = NCPatternSearch(x0, obj_fun, nlcon_fun, [], A, l, u,...
                            options.sl_ini, options.cvg_par, options.se_par, options.pattern);
    else
        error('This options is not finshed yet.');
    end

end