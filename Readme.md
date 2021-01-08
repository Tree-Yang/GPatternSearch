# A Family of MATLAB  Functions For Solving Optimization Problem Using Pattern Search



**ATTENTION: THIS PROJECT IS CREATED FOR STUDY PURPOSE. THEREFORE, THE CODE IS NOT THOROUGHLY TESTED!**

## Unconstrained Optimization

### Function name

```UCPatternSearch```

### Calling method

```[x_hist, f_hist] = UCPatternSearch(x0, obj_fun, sl_ini, cvg_par, se_par, pattern)```

where ```x_hist``` and ```f_hist```are the iteration histories of decision vector and objective function, respectively; ```x0```is initial decision vector; ```sl_ini```is initial step length for search; ```cvg_par```is a vector of parameters as convergent criteria;  ```se_par```is a vector controlling the shrinking and extending of search step length; and ```pattern```is a string indicating which pattern to use.

```'HJ'``` and ```'CS'```are possible value for ```pattern```, where ```HJ```means the pattern used by Hooke and Jeeves (1961), and ```'CS'``` refers to the 'Coordinate Search' pattern (Torczon, 1997).

Please refer to script "pstester.m" for proposed parameter values and examples of calling the function.

### References

Torczon V. (1997) On the convergence of pattern search, SIAM J. Optim, 7(1): 1-25.

Hooke R., Jeeves T.A. (1961) "Direct search" solution of numerical and statistical problems, J. Assoc. Comput. Mach, 8: 212-229.

## Bound Constrained Optimization

### Function name

``` PatternSearch```

### Calling method

```[x_hist, f_hist] = BCPatternSearch(x0, obj_fun, bounds ,sl_ini, cvg_par, se_par, pattern)```

where the inputs and outputs have the same meaning with function ```UCPatternSearch```, and ```bounds```contains the lower and upper bounds of decision variables:

$\mathbf{lb} \le \mathbf x \le \mathbf{ub}$

where ```lb = bounds(:,1)```and ```ub = bounds(:,2)```.

Please refer to script "pstester.m" for proposed parameter values and examples of calling the function.

 ### References

Lewis R.M., Torczon V. (1999) Pattern search algorithms for bound constrained minimization, SIAM J. Optim., 9(4): 1082-1099.

## Linearly Constrained Optimization

### Function name

```LCPatternSeach```

### Calling method

```[x_hist, f_hist] = LCPatternSearch(x0, obj_fun, A, l, u, sl_ini, cvg_par, se_par)```

where the inputs and outputs have the same meaning with function ```UCPatternSearch```and ```BCPatternSearch```,  and ```A```, ```l```, and ```u```are used to define linear constraint:

$\mathbf{l} \le \mathbf {Ax} \le \mathbf{u}$

Please refer to script "pstester.m" for proposed parameter values and examples of calling the function.

### References

Lewis R.M., Torczon V. (2000) Pattern Search methods for linearly constrained minimization, SIAM J. Optim., 10(3): 917-941.

## Nonlinearly Constrained Optimization

### Function name

```NCPatternSearch```

### Calling method

```[x_hist, f_hist, c_hist] = NCPatternSearch(x0, obj_fun, nlcon_fun, bounds, A, l, u, sl_ini, cvg_par, se_par, pattern)```

where the inputs and outputs have the same meaning with function ```UCPatternSearch```, ```BCPatternSearch``` and ```LCPatternSearch```;  ```c_histand```is the iteration history of nonlinear constraints, and  ```nlcon_fun``` is the function handle of nonlinear constraints.

Please refer to script "pstester.m" for proposed parameter values and examples of calling the function.

**ATTENTION: When nonlinear constraint are involved, the 'CS' pattern usually work more stable than the 'HJ' pattern.**

### References

Lewis R.M., Torczon V. (2002) A globally convergent augmented Lagrangian pattern search algorithm for optimization with general constraints and simple bounds, SIAM J. Optim, 12(4): 1075-1089.

Conn A.R., Gould N.I.M, Toint P.L. (1991) A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds, SIAM Numer. Anal., 28(2): 545-572.



## A Unified Solver

### Function name

```PatternSearchSolver```

### Calling method

```[x_hist, f_hist, c_hist] = PatternSearchSolver(x0, obj_fun, bounds, A, l, u, nlcon_fun,options)```

where ```options```is a structure containing algorithm settings, whose fields include ```sl_ini```, ```cvg_par```, ```se_par```, and ```pattern```.

Please refer to script "pstester.m" for proposed parameter values and examples of calling the function.





