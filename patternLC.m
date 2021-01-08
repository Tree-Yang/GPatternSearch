%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: define the pattern
%-----------------------------------------------
function patternMat = patternLC(A, l, u, xk, eps0)
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