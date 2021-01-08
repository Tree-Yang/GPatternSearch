
%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: bi-section
%-----------------------------------------------
function [x_fea,f_fea,nlc_fea] = bisec(ite, x1, f1, nlc1, x0, f0, nlc0, ...
                                    obj_fun, nlcon_fun, bisecn_max, bisec_eps)
    %Yang, JS; 2020-08-11

    %INPUT LIST
    %----------------------------------------
    %ite       : number of iteration (sub-problem)
    %x1        : a infeasible point
    %f1        : objective functionvalye at x1
    %nlc1      : value of nonlinear constraint functions at x1
    %x0        : a feasible point
    %f0        : objective function value at x0
    %nlc0      : value of nonlinear constraint functions at x0
    %obj_fun   : function handle of objective function
    %nlcon_fun : function handle of nonlinear constraint function
    %bisecn_max: maximum number of bisection
    %bisec_eps : exit tolerence of bisection

    %OUTPUT LIST
    %----------------------------------------
    %x_fea  : a feasible point by bisection
    %f_fea  : objective function value at x_fea
    %nlc_fea: nonlinear constraint function value at x_fea
    
    global nTrialGPS

    %number of variable
    n_x       = length(x0);
    %number of nonlinear constraints
    n_nlc     = length(nlc0);

    %stack of points and function values
    x_bschist        = zeros(n_x, bisecn_max+1);
    f_bschist        = zeros(1, bisecn_max+1);
    nlc_bschist      = zeros(n_nlc, bisecn_max+1);
    x_bschist(:,1)   = x1;
    f_bschist(:,1)   = f1;
    nlc_bschist(:,1) = nlc1;


    x_l   = x0;
    f_l   = f0;
    nlc_l = nlc0;
    x_r   = x1;
    f_r   = f1;
    nlc_r = nlc1;

    for kb = 1:1:bisecn_max

        fprintf('------------------------------------------------------------->\n');
        fprintf('Bisectioin Number:\t %d\n', kb);
        fprintf('------------------------------------------------------------->\n');

        %print every bisection point
        fprintf('Current point of bisection:\n');
        valuedisplay(x_r, 'x', 5);
        fprintf('Current Objective of bisection:\n');
        fprintf('f_o = %15.6f\n', f_r);
        fprintf('Current nonlinear constraints of bisection:\n');
        valuedisplay(nlc_r, 'nlc', 5);

        if norm(x_r - x_l) < bisec_eps
            fprintf(['The bisection procedure break since ',...
                        'The start point is too close to the end point!\n']);
            x_fea               = x_l;
            f_fea               = f_l;
            nlc_fea             = nlc_l;
            x_bschist(:,kb+1)   = x_l;
            f_bschist(:,kb+1)   = f_l;
            nlc_bschist(:,kb+1) = nlc_l;
            break;
        end

        %bisection
        x_c                 = (x_l+x_r)/2;
        f_c                 = obj_fun(x_c);
        nlc_c               = nlcon_fun(x_c);
        x_bschist(:,kb+1)   = x_c;
        f_bschist(:,kb+1)   = f_c;
        nlc_bschist(:,kb+1) = nlc_c;
        nTrialGPS           = nTrialGPS + 1;

        if sum(nlc_c > 0) == 0 && f_c < f_l
            fprintf(['The bisection procedure break since ',...
                        'a descent feasible point is found!\n']);
            x_fea   = x_c;
            f_fea   = f_c;
            nlc_fea = nlc_c;
            break;
        else
            x_r = x_c;
            % f_r = f_c;
            % nlc_r = nlc_c;
        end
    end

    if kb == bisecn_max
        fprintf(['The bisection procedure break since ',...
                    'too many bisection has been made!\n']);
        x_fea               = x_l;
        f_fea               = f_l;
        nlc_fea             = nlc_l;
        x_bschist(:,kb+1)   = x_l;
        f_bschist(:,kb+1)   = f_l;
        nlc_bschist(:,kb+1) = nlc_l;
    end

    %print the final bisection point
    fprintf('Final point of bisection:\n');
    valuedisplay(x_fea, 'x', 5);
    fprintf('Final Objective of bisection:\n');
    fprintf('f_o = %15.6f\n', f_fea);
    fprintf('Final nonlinear constraints of bisection:\n');
    valuedisplay(nlc_fea, 'nlc', 5);

    %save the iteration history of bisection into files
    x_bschist   = x_bschist(:,kb+1);
    f_bschist   = f_bschist(:,kb+1);
    nlc_bschist = nlc_bschist(:,kb+1);
    bschistfn   = ['SubOptimization_',num2str(ite),'_Bisection.mat'];
    save(bschistfn, 'x_bschist', 'f_bschist', 'nlc_bschist', '-v7.3');

end