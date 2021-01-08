%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: poll for pattern directions
%-----------------------------------------------
function [x, f, nlc, alf, flag] = polls(x0, f0, nlc0, sl, obj_fun, lb, ub, ...
                                        nlcon_fun, S, lambda_k, mu_k)
    %Yang, JS; 2020-08-09

    global nTrialGPS
    
    pollOpt   = 'full';
    % pollOpt   = 'first';
    
    %poll flag: 1 for successful poll; 0 for unsuccessful poll
    flag      = 0;
    %number of variables
    n_x       = length(x0);
    %matrix of pattern direction (positive,half)
    patternMat= eye(n_x);

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
    
    for ii = 1:1:n_x
        %trial points
        x_trial       = x - sl*patternMat(:,ii);
        if sum(x_trial < lb)==0 && sum(x_trial > ub)==0
            %function value at trial points
            f_trial   = obj_fun(x_trial);
            nlc_trial = nlcon_fun(x_trial);
            alf_trial = augLangF(f_trial, nlc_trial);
            ctr_trial = ctr_trial+1;
            nTrialGPS = nTrialGPS + 1;
        end

        if alf_trial >= alf
            %trial points
            x_trial       = x + sl*patternMat(:,ii);
            if sum(x_trial < lb)==0 && sum(x_trial > ub)==0 
                %function value at trial points
                f_trial   = obj_fun(x_trial);
                nlc_trial = nlcon_fun(x_trial);
                alf_trial = augLangF(f_trial, nlc_trial);
                ctr_trial = ctr_trial+1;
                nTrialGPS = nTrialGPS + 1;
            end
        end
        
        %update decision variables and function values when the polling succeed
        if alf_trial < alf && sum(x_trial < lb)==0 && sum(x_trial > ub)==0
            x   = x_trial;
            f   = f_trial;
            nlc = nlc_trial;
            alf = alf_trial;
        end
        
        if ~strcmpi(pollOpt, 'full')
            break;
        end

    end

    if alf < alf0
        flag = 1;
    end

    %output the number of trial
    fprintf('------------------------------------------------------------->\n')
    fprintf('Number of trial in the current step is %d\n', ctr_trial);
    fprintf('------------------------------------------------------------->\n')
    
end