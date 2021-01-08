%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: Linearly constrained pattern search for sub-problems
%-----------------------------------------------
function [x_hist, f_hist, nlc_hist, alf_hist] = LCPatternSearchSub(x0, f0, nlc0, obj_fun, ...
                                    A, l, u, nlcon_fun, sl_ini, S, lambda_k,...
                                    mu_k ,cvg_par, se_par)
    %Yang Jiashu; 2021-01-08
    
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');
    fprintf('Solution of the sub-problem:\n');
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');

    %number of variables
    n_x     = length(x0); 
    %number of nonlinear constraint function
    n_nlc    = length(nlc0);

    %augmented lagrangian function
    sii      = diag(S);
    siii     = (1./sii)';
    siid     = sii/mu_k;
    mu_siii  = mu_k/2*siii;
    zeronlc  = zeros(n_nlc,1);
    augLangF = @(of, nlcf) of + mu_siii*((max(zeronlc,...
                            lambda_k+siid.*nlcf)).^2-lambda_k.^2);
        
    %convergent parameters
    sl_esp   = cvg_par(1);
    ite_max  = cvg_par(2);

    %initialization
    sl       = sl_ini;
    ite      = 0;
    x        = x0;
    f        = f0;
    nlc      = nlc0;
    alf      = augLangF(f, nlc);

    x_hist   = zeros(n_x, ite_max);
    f_hist   = zeros(1, ite_max);
    nlc_hist = zeros(1, ite_max);
    alf_hist = zeros(1, ite_max);

    fprintf('Initial design:\n');
    valuedisplay(x, 'x', 5)
    fprintf('Initial function value:\n');
    fprintf('f_val = %15.6f\n', f);
    fprintf('Initial nonlinear constraint value:\n');
    valuedisplay(nlc, 'nlc', 5)
    fprintf('Initial value of augmented Lagrangian function:\n');
    fprintf('alf_val = %15.6f\n', alf);

    x_hist(:,1)   = x;
    f_hist(:,1)   = f;
    nlc_hist(:,1) = nlc;
    alf_hist(:,1) = alf;

    %start iteration
    while sl >= sl_esp && ite <= ite_max

        ite = ite + 1;

        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');
        fprintf('Pattern search number: %d\n', ite);
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');

        eps0 = 1e-3;
        patternMat   = patternLC(A, l, u, x, eps0);

        [x, f, nlc, alf, flag] = pollsLC(x, f, nlc, sl, obj_fun, A, l, u,...
                                    patternMat, nlcon_fun, S, lambda_k, mu_k);

        %record the point and function value
        x_hist(:,ite+1)   = x;
        f_hist(:,ite+1)   = f;
        nlc_hist(:,ite+1) = nlc;
        alf_hist(:,ite+1) = alf;
        fprintf('Current design:\n');
        valuedisplay(x, 'x', 5);
        fprintf('Current function value:\n');
        fprintf('f_val = %15.6f\n', f);
        fprintf('Current nonlinear constraint value:\n');
        valuedisplay(nlc, 'nlc', 5)
        fprintf('Current value of augmented Lagrangian function:\n');
        fprintf('alf_val = %15.6f\n', alf);

        %update the parameters of step length
        sl = slUpdate(sl, flag, se_par);
    end

    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');
    fprintf('The step length parameter is %15.6f\n', sl);
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++>\n');
    if ite >= ite_max
        fprintf(['Pattern Search terminated since the iteration number ', ...
                    'reached the maximum iteration number!\n']);
    else
        fprintf(['Pattern Search terminated since the step length parameter ',...
                    'is less than the tolerence!\n']);
    end
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++<\n');

    x_hist   = x_hist(:,1:ite+1);
    f_hist   = f_hist(:,1:ite+1);
    nlc_hist = nlc_hist(:,1:ite+1);
    alf_hist = alf_hist(:,1:ite+1);

end