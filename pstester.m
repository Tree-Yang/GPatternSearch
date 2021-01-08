clc; clear;
fun     = @psobj;
nonlcon = @ellipsetilt;
options = optimoptions('patternsearch','Display','iter');
sl_ini  = 1.0;
cvg_par = [1e-5, 100];
se_par  = [0.5,1.2];

%===================================================================================================
%Unconstrained pattern search
%------------------------------------------------------------------------------
% Direct solve
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% [x_hist, f_hist] = UCPatternSearch(x0, fun, sl_ini, cvg_par, se_par, 'CS');

% Call solver
%+++++++++++++++++++++++++++++++++++
% options1.sl_ini   = sl_ini;
% options1.cvg_par  = cvg_par;
% options1.se_par   = se_par;
% options1.pattern  = 'HJ';
% x0 = [-0;-0];
% [x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,[],[],[],[],[],options1);

% Bulit-in
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% x = patternsearch(fun,x0,[],[],[],[],[],[],[],options);
%===================================================================================================

%===================================================================================================
%Bound constrained pattern search
%------------------------------------------------------------------------------
% Direct solve
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% bounds = [-0.5, 0.5; -1, 1];
% [x_hist, f_hist] = BCPatternSearch(x0, fun, bounds, sl_ini, cvg_par, se_par, 'CS');

% Call solver
%+++++++++++++++++++++++++++++++++++
% options2.sl_ini   = sl_ini;
% options2.cvg_par  = cvg_par;
% options2.se_par   = se_par;
% options2.pattern  = 'HJ';
% bounds            = [-0.5, 0.5; -1, 1];
% x0 = [-0;-0];
% [x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,bounds,[],[],[],[],options2);

% Bulit-in
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% lb = [-0.5;-1];
% ub = [0.5;1];
% x = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);
%===================================================================================================

%===================================================================================================
%Linear constrained pattern search
%------------------------------------------------------------------------------
% Direct solve
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% A = [3,2; 4,7];
% u = [-3;1];
% l = -[100; 100];
% [x_hist, f_hist] = LCPatternSearch(x0, fun, A, l, u, sl_ini, cvg_par, se_par);

% Call solver
%+++++++++++++++++++++++++++++++++++
% options3.sl_ini   = sl_ini;
% options3.cvg_par  = cvg_par;
% options3.se_par   = se_par;
% options3.pattern  = 'HJ';
% x0 = [-0;-0];
% A = [3,2; 4,7];
% u = [-3;1];
% l = -[100; 100];
% [x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,[],A,l,u,[],options3);

% Bulit-in
%+++++++++++++++++++++++++++++++++++
% x0 = [-0;-0];
% A = [3,2; 4,7];
% b = [-3;1];
% x = patternsearch(fun,x0,A,b,[],[],[],[],[],options);
%===================================================================================================

%===================================================================================================
%Nonlinear constrained pattern search I
%------------------------------------------------------------------------------
% Direct solve I
%+++++++++++++++++++++++++++++++++++
% bounds = [-10,10;-10,10];
% x0     = [-1;0.5];
% [x_hist, f_hist, c_hist] = NCPatternSearch(x0, fun, nonlcon, bounds,...
%                                                 [], [], [], sl_ini, cvg_par, se_par, 'CS');

% Direct solve II
%+++++++++++++++++++++++++++++++++++
% x0     = [-1;0.5];
% A      = eye(2);
% l      = [-10; -10];
% u      = [10; 10];
% [x_hist, f_hist, c_hist] = NCPatternSearch(x0, fun, nonlcon, [],...
%                                                 A, l, u, sl_ini, cvg_par, se_par, 'CS');

% Call solver I
%+++++++++++++++++++++++++++++++++++
% options4.sl_ini   = sl_ini;
% options4.cvg_par  = cvg_par;
% options4.se_par   = se_par;
% options4.pattern  = 'CS';
% bounds = [-10,10;-10,10];
% x0     = [-1;0.5];
% [x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,bounds,[],[],[],nonlcon,options4);

% Call solver II
%+++++++++++++++++++++++++++++++++++
options5.sl_ini   = sl_ini;
options5.cvg_par  = cvg_par;
options5.se_par   = se_par;
options5.pattern  = 'CS';
x0     = [-1;0.5];
A      = eye(2);
l      = [-10; -10];
u      = [10; 10];
[x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,[],A,l,u,nonlcon,options5);

% Bulit-in
%+++++++++++++++++++++++++++++++++++
% x0     = [-1;0.5];
% lb     = [-10;-10];
% ub     = [10; 10];
% x = patternsearch(fun,x0,[],[],[],[],lb,ub,nonlcon,options);
%===================================================================================================

%===================================================================================================
%Nonlinear constrained pattern search I
%------------------------------------------------------------------------------
% Direct solve II
%+++++++++++++++++++++++++++++++++++
% x0= [-1;0.5];
% A = [3,2; 4,7];
% u = [-3;1];
% l = -[100; 100];
% [x_hist, f_hist, c_hist] = NCPatternSearch(x0, fun, nonlcon, [],...
%                                                 A, l, u, sl_ini, cvg_par, se_par, 'CS');

% Call solver II
%+++++++++++++++++++++++++++++++++++
% options6.sl_ini   = sl_ini;
% options6.cvg_par  = cvg_par;
% options6.se_par   = se_par;
% options6.pattern  = 'CS';
% x0= [-1;0.5];
% A = [3,2; 4,7];
% u = [-3;1];
% l = -[100; 100];
% [x_hist, f_hist, c_hist] = PatternSearchSolver(x0,fun,[],A,l,u,nonlcon,options6);

% Bulit-in
%+++++++++++++++++++++++++++++++++++
% x0= [-1;0.5];
% A = [3,2; 4,7];
% b = [-3;1];
% x = patternsearch(fun,x0,A,b,[],[],[],[],nonlcon,options);
%===================================================================================================

% save('WorkSpace.mat','-v7.3');

%===================================================================================================
function y = psobj(x)
    y = exp(-x(1)^2-x(2)^2)*(1+5*x(1) + 6*x(2) + 12*x(1)*cos(x(2)));
end

function [c,ceq] = ellipsetilt(x)
    ceq = [];
    c = x(1)*x(2)/2 + (x(1)+2)^2 + (x(2)-2)^2/2 - 2;
end