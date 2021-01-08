%==========================================================================
%SUBFUNCTION FOR PATTERN SEARCH: update of step length parameter
%-----------------------------------------------
function sl = slUpdate(sl, flag, se_par)

    %Yang, JS; 2020-08-08

    shrink = se_par(1);
    extend = se_par(2);
    
    if shrink >= 1
        error('The shrink parameter of step length should be smaller than 1.0!');
    end
    
    if extend < 1
        error('The shrink parameter of step length should not be smaller than 1.0!');
    end

    if flag == 0
        % when the polling fails
        sl = sl*shrink;
    else
        % when the polling succeeds
        sl = sl*extend;
    end

end
%==========================================================================