%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: poll for pattern directions
%-----------------------------------------------
function [x, f, nlc, alf, flag] = pollsLC(x0, f0, nlc0, sl, obj_fun, A, l, u,...
                                    patternMat, nlcon_fun, S, lambda_k, mu_k)
    %Yang, JS; 2020-08-08
    
    global nTrialGPS
    
    %poll flag: 1 for successful poll; 0 for unsuccessful poll
    flag       = 0;

    %number of pattern direction
    n_d        = size(patternMat,2);

    %number of nonlinear constraint function
    n_nlc     = length(nlc0);

    %augmented lagrangian function
    sii       = diag(S);
    siii      = (1./sii)';
    siid      = sii/mu_k;
    mu_siii   = mu_k/2*siii;
    zeronlc   = zeros(n_nlc,1);
    augLangF  = @(of, nlcf) of + mu_siii*((max(zeronlc,...
                            lambda_k+siid.*nlcf)).^2-lambda_k.^2);
    
    %trial and polls
    x         = x0;
    f         = f0;
    nlc       = nlc0;
    alf0      = augLangF(f, nlc);
    alf       = alf0;
    ctr_trial = 0;

    f_trial   = f;
    nlc_trial = nlc;
    alf_trial = alf;
    
    for ii = 1:1:n_d
        %trial points
        x_trial   = x0 + sl*patternMat(:,ii);
        if sum(A*x_trial > u) == 0 && sum(A*x_trial < l) ==0
            %function value at trial points
            f_trial   = obj_fun(x_trial);
            nlc_trial = nlcon_fun(x_trial);
            alf_trial = augLangF(f_trial, nlc_trial);
            ctr_trial = ctr_trial+1;
            nTrialGPS = nTrialGPS + 1;
        end
        
        if alf_trial < alf
            x   = x_trial;
            f   = f_trial;
            nlc = nlc_trial;
            alf = alf_trial;
            flag = 1;
            break;
        end

    end

    %output the number of trial
    fprintf('------------------------------------------------------------->\n')
    fprintf('Number of trial in the current step is %d\n', ctr_trial);
    fprintf('------------------------------------------------------------->\n')
    
end