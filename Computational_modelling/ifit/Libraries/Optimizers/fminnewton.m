function [pars,fval,exitflag,output] = fminnewton(varargin)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = fminnewton(FUN,PARS,[OPTIONS],[CONSTRAINTS], ...) Newton search
%
% This minimization method uses the Dogleg trust region, Newton model, dense algorithm 
% The objective function has syntax: criteria = objective(p)
% 
% Calling:
%   fminnewton(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminnewton(fun, pars, options) same as above, with customized options (optimset)
%   fminnewton(fun, pars, options, fixed) 
%     is used to fix some of the parameters. The 'fixed' vector is then 0 for
%     free parameters, and 1 otherwise.
%   fminnewton(fun, pars, options, lb, ub) 
%     is used to set the minimal and maximal parameter bounds, as vectors.
%   fminnewton(fun, pars, options, constraints) 
%     where constraints is a structure (see below).
%   fminnewton(problem) where problem is a structure with fields
%     problem.objective:   function to minimize
%     problem.x0:          starting parameter values
%     problem.options:     optimizer options (see below)
%     problem.constraints: optimization constraints
%   fminnewton(..., args, ...)
%     sends additional arguments to the objective function
%       criteria = FUN(pars, args, ...)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminnewton(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string): criteria = FUN(PARS)
%  It needs to return a single value or vector.
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess. PARS can also be given as a single-level structure.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%     o=fminnewton('defaults')
%  options.MinFunEvals sets the minimum number of function evaluations to reach
%  An empty OPTIONS sets the default configuration.
%
%  CONSTRAINTS may be specified as a structure
%   constraints.min=   vector of minimal values for parameters
%   constraints.max=   vector of maximal values for parameters
%   constraints.fixed= vector having 0 where parameters are free, 1 otherwise
%   constraints.step=  vector of maximal parameter changes per iteration
%   constraints.eval=  expression making use of 'p', 'constraints', and 'options' 
%                        and returning modified 'p'
%                      or function handle p=@constraints.eval(p)
%  An empty CONSTRAINTS sets no constraints.
%
%  Additional arguments are sent to the objective function.
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
%
% Reference: W. Press, Numerical Recipes, Cambridge (1988)
% Contrib: C. T. Kelley, 1998, Iterative Methods for Optimization [ntrust]
%
% Version: Aug. 22, 2017
% See also: fminsearch, optimset
% (c) E.Farhi, ILL. License: EUPL.

% default options for optimset
if nargin == 0 || (nargin == 1 && strcmp(varargin{1},'defaults'))
  options=optimset; % empty structure
  options.Display='';
  options.TolFun =1e-3;
  options.TolX   =1e-8;
  options.MaxIter='10*numberOfVariables';
  options.MaxFunEvals=1000;
  options.algorithm  = [ 'Steihaug Newton-CG-Trust (by Kelley) [' mfilename ']' ];
  options.optimizer = mfilename;
  pars = options;
  return
end

[pars,fval,exitflag,output] = fmin_private_wrapper(mfilename, varargin{:});

